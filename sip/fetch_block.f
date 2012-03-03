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
      subroutine fetch_block(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, direct_flag)
c--------------------------------------------------------------------------
c   Fetches the current block of the array listed in the result_array
c   field of the operation.
c
c   The array must be declared as a distributed array.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      logical direct_flag

      integer i, j, k
      integer array, array_type, proc, ierr
      integer blk, blkndx, request, old_request
      integer size
      integer find_block_destination
      integer get_block_request
      
      integer pst_get_company_comm
      integer comm
      integer flag

      array = op(c_result_array)
      array_type = array_table(c_array_type,array)
      if (array_type .ne. distributed_array) then
         print *,'Error in fetch_block: Array ',array,' should be ',
     *           'declared as a distributed array.'
         call abort_job()
      endif
      
      proc = find_block_destination(array, array_table,
     *               narray_table, index_table, nindex_table,
     *               block_map_table, nblock_map_table)
      if (proc .lt. 0) then
         print *,'Error: Current block of array ',array,' was not',
     *            ' found in block_map_table.'
         print *,'Current indices: ',(index_table(c_current_seg,i),
     *             i = 1, nindex_table)
         print *,'Dumping block map table:'
         do i = 1, nblock_map_table
            print *,'Entry ',i,': ',(block_map_table(j,i),
     *               j=1,lblock_map_entry)
         enddo
         call dump_block_ids()
         call abort_job()
      endif

      call create_current_block(array,array_table, narray_table,
     *                 index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., direct_flag, blk, ierr)
      if (ierr .le. 0) then
         print *,'Error: Cannot create block for array ',array
         print *,'ierr = ',ierr,' blk = ',blk
         call dump_block_ids()
         call abort_job()
      else
         blkndx = ierr
      endif 

c      call check_block_consistency(array, blk, blkndx)

      call get_block_computed_flag(array, blk, blkndx, flag)
      if (flag .eq. 0) then
         call get_block_created_flag(array, blk, blkndx, flag)
         if (flag .eq. 0) then

c----------------------------------------------------------------------------
c   Turn on the block_computed_flag.  This tells us the block was put into
c   use on this instruction.
c----------------------------------------------------------------------------

            call set_opblock(array, blk, blkndx, op)
            call set_block_computed_flag(array, blk, blkndx, 1)

c----------------------------------------------------------------------------
c   If this block was the output of a previous loop, the scrub flag will be 
c   set.  If this is the case, we must turn it off here.  That prevents it
c   from being scrubbed again and used for something else.
c----------------------------------------------------------------------------

            call set_block_scrub_flag(array, blk, blkndx, 0)
         endif
      endif

c         call get_blk_header(flag, blkndx, c_block_flags)
c         print 10101,me,array,current_line,blkndx,flag,proc
c10101 format('Task ',i2,' FETCH array ',i5,' line ',i5,
c     *         ' blkndx ',i5,
c     *          ' flag ',i5,
c     *          ' proc ',i2))

      if (proc .eq. my_company_rank) then
c         call peg_block_access_stats(c_local_read)
         return   ! block resides on this proc.
      endif

c-------------------------------------------------------------------------
c   Communication may be required.  First check for a previous communication 
c   that may still be in progress.
c--------------------------------------------------------------------------

      comm = pst_get_company_comm(me) 
      old_request = get_block_request(array,blk, blkndx)

      if (old_request .eq. mpi_request_null) then

c--------------------------------------------------------------------------
c   Request the block from its host processor.
c   The request is only made if the block has not been previously requested.
c   A "wait" must be issued elsewhere before the block may be used.
c--------------------------------------------------------------------------


c---------------------------------------------------------------------------
c   Check the block_persistence_flag.  If it is already turned on, the
c   block is available from a previous access, and there is no need to
c   request it again.
c---------------------------------------------------------------------------

            call get_block_persistence_flag(array, blk, blkndx, flag)
            if (flag .eq. 0) then
               call get_actual_blocksize(array, blk, blkndx,
     *              array_table, narray_table,
     *              index_table, nindex_table,
     *              segment_table, nsegment_table, size)

               call frequestblk(array, proc, blk, blkndx, size,
     *                      request, comm)

               call set_block_request(array, blk, blkndx, request) 
               if (request .ne. mpi_request_null) 
     *            call blkmgr_insert_block_in_list(comm_list_head,
     *                    comm_list_tail, blkndx, c_comm_list_ptr,
     *                       .true.)

c--------------------------------------------------------------------------
c   Set the block_persistence_flag.  This guarantees that the block will
c   remain available as long as possible before having to do an actual
c   new request from it's "home" processor.
c
c   Setting the flag also updates its age, which is used in the allocation
c   scheme.  
c--------------------------------------------------------------------------

               call set_block_persistence_flag(array, blk, blkndx, 1)
               call blkmgr_insert_block_in_list(persistent_list_head,
     *                persistent_list_tail, blkndx, c_persistent_ptr,
     *                .true.)
            endif
      endif

c         call get_blk_header(flag, blkndx, c_block_flags)
c         print 10102,me,array,current_line,blkndx,flag
c10102 format('Task ',i2,' FETCH RET array ',i5,' line ',i5,
c     *          ' blkndx ',i5,
c     *          ' flag ',i5,
c     *          ' proc ',i2))

      return
      end
