C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine stack_distribution(optable, noptable,
     *                    array_table, narray_table,
     *                    index_table, nindex_table,
     *                    segment_table, nsegment_table,
     *                    scalar_table, nscalar_table,
     *                    block_map_table, nblock_map_table,
     *                    proctab, scr, nscr, nblk_overhead,
     *                    stack_blocksizes, nstacks, stack_blocks,
     *                    stack_algorithm_type, ncompany_procs,
     *                    niocompany, ierr, error_print_flag)
      implicit none
      include 'interpreter.h'
      include 'machine_types.h'
      include 'blkmgr.h'
      include 'mpif.h' 
      include 'dbugcom.h'

      integer noptable, narray_table, nindex_table, nsegment_table,
     *        nscalar_table, nblock_map_table, nscr, nstacks
      integer nblk_overhead
      integer nbytes_per_blk, nbytes_per_group
      integer optable(loptable_entry,noptable)
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      double precision scalar_table(nscalar_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer proctab(2,*), scr(nstacks,*)
      integer stack_blocksizes(nstacks)
      integer stack_blocks(nstacks)
      integer stack_line(nstacks)
      integer stack_algorithm_type
      integer niocompany
      integer ierr
      logical error_print_flag

      integer me, nprocs, mpierr
      integer i, j,k, ncompany_procs, my_company_rank, comm
      integer pst_get_company_comm
      integer nleft, total_blocks
      integer*8 noverhead_bytes, i8temp
      integer noverhead_words
      integer ngroups, nwords_per_group, nblocks_per_stack
      integer maxblk, need, have, surplus

      integer nresidual, nresidual_blocks
      integer nresidual_old
      integer*8 nwadd, nbadd
      integer new_workers, nworkers
      integer*8 total_for_job, total_words 
      double precision weight

      ierr = 0
      call mpi_comm_rank(mpi_comm_world, me, mpierr)
      call mpi_comm_size(mpi_comm_world, nprocs, mpierr)
      if (dbg) print *,'Determine stack distribution for task ',me 
      comm =  pst_get_company_comm(me)
      call mpi_comm_rank(comm, my_company_rank, mpierr)

c--------------------------------------------------------------------------
c   Zero out the block counters.
c--------------------------------------------------------------------------

      if (nscr .lt. noptable*nstacks) then
          print *,'Error: insufficient scratch memory in ',
     *        'stack_distribution: nscr ',nscr,' need ',
     *         noptable*nstacks
          call abort_job()
      endif

      do j = 1, noptable
         do i = 1, nstacks
            scr(i,j) = 0
         enddo
      enddo

c---------------------------------------------------------------------------
c   Make a pass simulating the optable execution to get an estimate of 
c   number of blocks of each type needed at each instruction line.
c---------------------------------------------------------------------------

      call optable_loop_sim(optable, noptable, array_table,
     *                   narray_table,
     *                   index_table, nindex_table,
     *                   segment_table, nsegment_table, block_map_table,
     *                   nblock_map_table,
     *                   scalar_table, nscalar_table, proctab, 
     *                   stack_blocksizes, scr, nstacks)

c---------------------------------------------------------------------------
c   Determine the max number of blocks needed for each stack type.
c---------------------------------------------------------------------------

      total_blocks = 0
      do i = 1, nstacks
         stack_blocks(i) = 0
         stack_line(i)  = 0
         do j = 1, noptable
            if (scr(i,j) .gt. stack_blocks(i)) then
               stack_blocks(i) = scr(i,j)
               stack_line(i)   = optable(c_lineno,j)
c               print *,'Task ',me,' New max. for stack ',i,' instr. ',
c     *           j,' opcode ',optable(c_opcode,j),' line ',
c     *           optable(c_lineno,j),' blocks = ',stack_blocks(i)
            endif
         enddo 

         total_blocks = total_blocks + stack_blocks(i)
      enddo
   
c------------------------------------------------------------------------
c   Make sure each stack has at least 1 block available.
c------------------------------------------------------------------------

      do i = 1, nstacks
         if (stack_blocks(i) .eq. 0) then
            stack_blocks(i) = 1
            total_blocks = total_blocks + 1
         endif
      enddo
 
c-------------------------------------------------------------------------
c   Need 3 scratch blocks on the largest stack.
c-------------------------------------------------------------------------

      stack_blocks(nstacks) = stack_blocks(nstacks) + 3
      total_blocks = total_blocks + 3
      if (dbg) print *,'Task ',me,
     *     ' Stack distribution after optable_loop_sim: ',
     * (stack_blocks(i),i=1,nstacks),' total_blocks ',total_blocks

c---------------------------------------------------------------------------
c   Final recheck: Are there too many blocks?
c--------------------------------------------------------------------------

c      if (total_blocks .gt. max_blkmgr_blocks) then 
c        if (error_print_flag) then
c         print *,'Error: The current job parameters require a total',
c     *           ' of ',total_blocks,' which exceeds the limit of ',
c     *           max_blkmgr_blocks
c         print *,'Either different segmentation of occupied/virtual/AO',
c     *           ' indices must be used, or a larger number of ',
c     *           'processors is required to run this job.'
c         print *,'Another alternative is to choose a different SIAL',
c     *           ' algorithm.'
c         print *,'2D stack: ',stack_blocks(1),
c     *      ' blocks because of  line ',stack_line(1)
c         print *,'OOOO/OOOV: ',stack_blocks(2),
c     *      ' blocks because of  line ',stack_line(2)
c         print *,'OOOV/OOVV: ',stack_blocks(3),
c     *      ' blocks because of  line ',stack_line(3)
c         print *,'OOVV: ',stack_blocks(4),
c     *      ' blocks because of  line ',stack_line(4)
c         print *,'OOVV/OVVV: ',stack_blocks(5),
c     *      ' blocks because of  line ',stack_line(5)
c         print *,'OVVV: ',stack_blocks(6),
c     *      ' blocks because of  line ',stack_line(6)
c         print *,'VVVV: ',stack_blocks(7),
c     *      ' blocks because of  line ',stack_line(7)
c        endif
c        ierr = 1
c        return
c      endif

      if (dbg) print *,'Task ',me,' total_blocks, max_blkmgr_blocks ',
     *     total_blocks, max_blkmgr_blocks

c---------------------------------------------------------------------------
c  Do the blocks take up too much memory?
c---------------------------------------------------------------------------

      total_words = 0
      do i = 1, nstacks

c---------------------------------------------------------------------------
c   When computing total_words, we must use a 64-bit temporary variable in
c   order to avoid integer overflow due to a possibly large number.
c---------------------------------------------------------------------------

         i8temp = stack_blocksizes(i) + lblock_id_data
         i8temp = i8temp * stack_blocks(i)
         total_words = total_words + i8temp
      enddo

c---------------------------------------------------------------------------
c  Add the overhead due to block management.
c---------------------------------------------------------------------------

      noverhead_bytes = nblk_overhead * intsize * 
     *                   total_blocks
      noverhead_words = (noverhead_bytes + bytes_per_double - 1) /
     *                   bytes_per_double
      total_words = total_words + noverhead_words
      if (total_words .gt. nscr) then
        if (error_print_flag) then
         print *,'Error: The current job configuration will require ',
     *         'the following minimum configuration of stacks to run:'
         print *,'2D stack: ',stack_blocks(1),' blocks of ',
     *            stack_blocksizes(1),' words because of line ',
     *            stack_line(1)
         print *,'OOOO/OOOV: ',stack_blocks(2),' blocks of ',
     *            stack_blocksizes(2),' words because of line ',
     *            stack_line(2)
         print *,'OOOV/OOVV: ',stack_blocks(3),' blocks of ',
     *            stack_blocksizes(3),' words because of line ',
     *            stack_line(3)
         print *,'OOVV: ',stack_blocks(4),' blocks of ',
     *            stack_blocksizes(4),' words because of line ',
     *            stack_line(4)
         print *,'OOVV/OVVV: ',stack_blocks(5),' blocks of ',
     *            stack_blocksizes(5),' words because of line ',
     *            stack_line(5)
         print *,'OVVV: ',stack_blocks(6),' blocks of ',
     *            stack_blocksizes(6),' words because of line ',
     *            stack_line(6)
         print *,'VVVV: ',stack_blocks(7),' blocks of ',
     *            stack_blocksizes(7),' words because of line ',
     *            stack_line(7)
         print *,'Total words requred: ',total_words,
     *           ' actual number available: ',nscr
         print *,'Memory shortfall: ',
     *         ((total_words-nscr)*8)/(1024*1024),
     *          ' Mbytes'
        endif

        ierr = 1
        return
      endif
      
      if (dbg) print *,'Task ',me,' total_words, nscr ',
     *                   total_words, nscr,' stack algorithm ',
     *                   stack_algorithm_type

      if (stack_algorithm_type .eq. 1) then

c--------------------------------------------------------------------------
c   Divide out any remaining blocks evenly among the different stacks.
c--------------------------------------------------------------------------

         nwords_per_group = 0
         do i = 1, nstacks
            nwords_per_group = nwords_per_group + stack_blocksizes(i) +
     *                      lblock_id_data
         enddo

         nbytes_per_group = nwords_per_group * bytes_per_double + 
     *                      nstacks * nblk_overhead * intsize
         nblocks_per_stack = (nscr-total_words) * bytes_per_double / 
     *                        nbytes_per_group

         if (dbg) print *,'Task ',me,' Adding ',nblocks_per_stack,
     *            ' to each stack'
         do i = 1, nstacks
            stack_blocks(i) = stack_blocks(i) + nblocks_per_stack
         enddo
      else if (stack_algorithm_type .eq. 2) then

c--------------------------------------------------------------------------
c   Memory is proportionally allocated among the non-zero stacks, from the
c   stack with the largest block_size down.
c   This is more useful for scf-type distributions.
c--------------------------------------------------------------------------

         nresidual = nscr - total_words
         nresidual_old = nresidual

         do i = nstacks, 1, -1
            if (stack_blocks(i) .ne. 0 .and. 
     *          nresidual .gt. 0) then
               weight =float(stack_blocks(i)) / total_blocks
               i8temp  = weight * nresidual_old * bytes_per_double
               nwadd = i8temp
               nbytes_per_blk = (stack_blocksizes(i) + lblock_id_data) 
     *            *bytes_per_double + nblk_overhead * intsize
               nbadd  = i8temp / nbytes_per_blk 
               noverhead_words = (nblk_overhead * intsize* nbadd +
     *                       bytes_per_double - 1) / bytes_per_double
               nresidual = nresidual - nbadd * (stack_blocksizes(i) +
     *                     lblock_id_data) - noverhead_words
               if (nresidual .gt. 0) then
                  if (dbg) print *,'Task ',me,' Adding ',nbadd,
     *                 ' blocks to stack ',
     *                  i 
                  stack_blocks(i) = stack_blocks(i) + nbadd
               endif
            endif
         enddo 
      else
         print *,'Error: Invalid algorithm_type parameter in ',
     *       'stack_distribution routine.'
         print *,'Value is ',stack_algorithm_type
         call abort_job()
      endif

      return
      end
