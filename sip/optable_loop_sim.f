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
      subroutine optable_loop_sim(optable, noptable, array_table, 
     *                   narray_table,
     *                   index_table, nindex_table, 
     *                   segment_table, nsegment_table, block_map_table,
     *                   nblock_map_table, 
     *                   scalar_table, nscalar_table, proctab,
     *                   stack_blksizes, block_count, nstacks)
c--------------------------------------------------------------------------
c   Simulates the execution of the optable, accumulating a count of 
c   the number of blocks needed on each stack for each instruction.
c--------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_barrier_data.h'
      include 'scratchpad.h'

      common /load_balance/load_balance
      logical load_balance

      integer noptable, narray_table, nindex_table, nsegment_table
      integer nblock_map_table, nscalar_table
      integer optable(loptable_entry,noptable)
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer proctab(2,*)
      double precision scalar_table(nscalar_table)
      integer nstacks
      integer stack_blksizes(nstacks)
      integer block_count(nstacks,noptable)
 
      integer ierr
      integer iop, iopsave, iblk
      integer index, result_array, result_type
      integer opcode
      integer i, j, k, n
      integer running_count(nstacks)

      integer level, type, lop, lopbegin, lopend
      integer max_level
      parameter (max_level = 100)
      integer current(nstacks, max_level)
      integer instr(max_level)
      integer matchind(mx_array_index)
      integer save_ind(mx_array_index,narray_table)

      integer ind(mx_array_index), seg(mx_array_index)
      integer stack, which_stack
      integer array, nind, numblks, iblock_map, proc
      integer val1, val2
      integer allocate_table(narray_table)
      integer instruction

      integer jseg, iseg(mx_array_index)
      integer iblock, nseg(mx_array_index), bseg(mx_array_index),
     *        eseg(mx_array_index), iwild(mx_array_index),
     *        maxrange(mx_array_index), maxrange_seg(mx_array_index)
      integer nblock, nwild
      integer start_op, end_op
      integer lookup, block_map_lookup
      logical debug
      logical countit
      logical partial_create

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, me, ierr)

      trace = .false.
      load_balance = .false.
      debug = .false.
      simulator = .true.

      do i = 1, narray_table
         allocate_table(i) = 0
      enddo

      do i = 1, nstacks
         running_count(i) = 0
      enddo
c      running_count(nstacks) = 3   ! scratch arrays for contraction routines.

c--------------------------------------------------------------------------
c   Clear the scratchpad areas.
c--------------------------------------------------------------------------

      do i = 1, mx_scratchpad
         scratchpad(i) = 0
      enddo

      do i = 1, mx_fscratchpad
         fscratchpad(i) = 0.
      enddo

c--------------------------------------------------------------------------
c   Save original array table indices.
c--------------------------------------------------------------------------

      do i = 1, narray_table
         do j = 1, mx_array_index
            save_ind(j,i) = array_table(c_index_array1+j-1,i)
         enddo
      enddo

c--------------------------------------------------------------------------
c   Main processing loop.
c--------------------------------------------------------------------------

         start_op = 0
         end_op   = noptable
         iop = 1 
 1000    continue
            opcode = optable(c_opcode, iop)
            iopsave = iop
            current_op = iop  ! save operation pointer for debugging  purposes.
            current_line = optable(c_lineno,iop)

            if (opcode .eq. call_op) then
               call handle_call(optable, noptable, proctab, debug,
     *                start_op, end_op, iop)
            else if (opcode .eq. return_op) then
               call handle_return(optable, noptable, debug,
     *                start_op, end_op, iop)
            else if (opcode .eq. go_to_op) then 
               call handle_goto(optable, noptable, debug,
     *                start_op, end_op, iop)
            else if (opcode .eq. do_op .or.
     *               opcode .eq. enddo_op) then
c               call doloop(optable, noptable, iop, index_table, 
c     *                  nindex_table, array_table, narray_table,
c     *                  block_map_table, segment_table,
c     *                  nsegment_table, 
c     *                  debug, .false.,
c     *                  start_op, end_op)
            else if (opcode .eq. pardo_op .or.
     *               opcode .eq. endpardo_op) then
c               if (load_balance) then
c                  call pardo_loadb(optable, noptable, iop, index_table,
c     *                  nindex_table, array_table, narray_table,
c     *                  block_map_table,
c     *                  comm, debug, .false.,
c     *                  start_op, end_op)
c               else
c                  call pardo_loop(optable, noptable, iop, index_table,
c     *                  nindex_table, array_table, narray_table,
c     *                  block_map_table,
c     *                  comm, debug, .false.,
c     *                  start_op, end_op)
c               endif
            else if (opcode .eq. exit_op) then
c               call handle_exit(optable, noptable, debug,
c     *                array_table, narray_table,
c     *                index_table, nindex_table, 
c     *                block_map_table,    
c     *                start_op, end_op, iop)
            else if (opcode .eq. cycle_op) then
c               call handle_cycle(optable, noptable, debug,
c     *                start_op, end_op, iop)
            endif

            if (iop .gt. noptable) go to 2000
            if (iopsave .ne. iop) then
               go to 1000
            endif

            if (iop .lt. start_op .or.
     *          iop .gt. end_op) go to 900 

            current_op = iop  ! save operation pointer for debugging  purposes.
            current_line = optable(c_lineno,iop)

            do i = 1, nstacks
               block_count(i,iop) = running_count(i)
            enddo

c---------------------------------------------------------------------------
c   Check for the following instructions: CREATE, DELETE, ALLOCATE,
c   or DEALLOCATE.
c---------------------------------------------------------------------------

            partial_create = .false.   

            if (optable(c_opcode,iop) .eq. create_op .or.
     *          optable(c_opcode,iop) .eq. delete_op) then

c---------------------------------------------------------------------------
c   Determine the number of blocks for this array on the current processor.
c---------------------------------------------------------------------------

               array = optable(c_result_array,iop)
               iblock_map = array_table(c_block_map,array)
               numblks    = array_table(c_numblks,array)
               nind       = array_table(c_nindex,array)

               do i = 1, nind
                  ind(i) = array_table(c_index_array1+i-1,array)
               enddo

c---------------------------------------------------------------------------
c   Is this a partial create/delete?
c---------------------------------------------------------------------------

               partial_create = .false.
               do k = 1, nind
                  if (optable(c_ind1+k-1,iop) .eq. 
     *                             wildcard_indicator) then
                     partial_create = .true.
                  endif    
               enddo

               if (partial_create) go to 500
                  
c---------------------------------------------------------------------------
c   Loop over the number of blocks.
c---------------------------------------------------------------------------

               do i = 1, numblks
                  proc = block_map_table(c_processor,iblock_map+i-1)
                  if (proc .eq. my_company_rank) then
               
c----------------------------------------------------------------------------
c   Determine blocksize.
c----------------------------------------------------------------------------

                     n = 1
                     do j = 1, nind
                        seg(j) =  block_map_table(c_block_map_seg+j-1,
     *                                          iblock_map+i-1)
                        index_table(c_current_seg,ind(j)) = seg(j)
                        call get_index_segment(ind(j), seg(j), 
     *                             segment_table, nsegment_table, 
     *                             index_table, nindex_table, 
     *                             val1, val2)
                        if (val1 .le. 0 .or. val2 .le. 0) go to 900
                        n = n * (val2 - val1 + 1)
                     enddo

c----------------------------------------------------------------------------
c   Determine the proper stack.
c----------------------------------------------------------------------------

                     stack = which_stack(stack_blksizes, nstacks, n) 

c---------------------------------------------------------------------------
c   Increment (or decrement) the running block count.
c---------------------------------------------------------------------------

                     if (optable(c_opcode,iop) .eq. create_op) then
                        block_count(stack,iop) = 
     *                            block_count(stack,iop)+1
                     else
                        block_count(stack,iop) =
     *                            block_count(stack,iop) - 1
                     endif
                  endif
               enddo 
            endif   ! create_op/delete_op
  500       continue
 
            if (optable(c_opcode,iop) .eq. allocate_op .or.
     *          (optable(c_opcode,iop) .eq. create_op .and. 
     *          partial_create)) then 

c----------------------------------------------------------------------------
c   Since we don't know the actual distribution of the non-wildcard indices
c   until runtime, we cannot directly calculate the size of each block.
c   Instead, we calculate the number of wildcard blocks, and use the maximum
c   size of each non-wildcard index to estimate an upper limit of the
c   size of each block.
c----------------------------------------------------------------------------

               array = optable(c_result_array,iop)
               allocate_table(array) = iop   ! save the instruction for dealloc

               nblock = 1
               nwild  = 0
               nind   = 0
               do i = 1, mx_array_index
                  if (array_table(c_index_array1+i-1,array) .gt. 0) then
                     nind = nind + 1
                     ind(i) = array_table(c_index_array1+i-1,array)
                     nseg(i) = index_table(c_nsegments, ind(i))
                     bseg(i) = index_table(c_bseg, ind(i))
                     eseg(i) = index_table(c_eseg, ind(i))

                     if (optable(c_ind1+i-1,iop) .eq. 
     *                       wildcard_indicator) then
                        nblock = nblock*nseg(i)
                        nwild = nwild + 1
                        iwild(nwild) = i
                        iseg(iwild(nwild)) = bseg(iwild(nwild)) 
                     endif
                  endif
               enddo

               do j = 1, nind
                  maxrange(j) = 0
                  maxrange_seg(j) = 0
                  do jseg = bseg(j), eseg(j)
                     call get_index_segment(ind(j),jseg,
     *                             segment_table, nsegment_table,
     *                             index_table, nindex_table,
     *                             val1, val2)
                     if (val1 .le. 0 .and. val2 .le. 0) go to 900
                     n = val2 - val1 + 1
                     if (n .gt. maxrange(j)) then 
                        maxrange(j) = n
                        maxrange_seg(j) = jseg
                     endif 
                  enddo
               enddo

               do iblock = 1, nblock
         
c--------------------------------------------------------------------------
c   Determine the blocksize.  First loop over the non-wildcard indices, 
c   using the maxrange for the size of each of their segments.
c--------------------------------------------------------------------------

                  n = 1
                  do i = 1, nind
                     if (optable(c_ind1+i-1,iop) .ne. 
     *                          wildcard_indicator) then
                        n = n * maxrange(i)
                        seg(i) = maxrange_seg(i)
                     endif
                  enddo

c---------------------------------------------------------------------------
c   Next, loop over the wildcard indices, determining the actual segment
c   size for the next block.
c---------------------------------------------------------------------------

                  do i = 1, nwild
                     call get_index_segment(ind(iwild(i)),
     *                                      iseg(iwild(i)),
     *                             segment_table, nsegment_table,
     *                             index_table, nindex_table,
     *                             val1, val2)
                     n = n * (val2 - val1 + 1) 
                     seg(i) = iseg(iwild(i))
                  enddo

c--------------------------------------------------------------------------
c   Determine the stack on which to place the block.
c--------------------------------------------------------------------------

                  stack = which_stack(stack_blksizes, nstacks, n)

c---------------------------------------------------------------------------
c   Adjust the block count for the stack.
c---------------------------------------------------------------------------

                  if (partial_create) then

c---------------------------------------------------------------------------
c   Only count it if the block_map_table indicates the block resides on
c   this processor.
c---------------------------------------------------------------------------

                     lookup =  block_map_lookup(seg, nind, array,
     *                       array_table(1,array),
     *                       index_table, nindex_table) 
                     proc = block_map_table(c_processor,lookup)
                     if (proc .eq. my_company_rank) 
     *                   block_count(stack,iop) =
     *                        block_count(stack,iop) + 1  
                  else
                     block_count(stack,iop) = 
     *                        block_count(stack,iop) + 1
                  endif

c--------------------------------------------------------------------------
c   Increment the "wildcard" segments.
c--------------------------------------------------------------------------

                  if (nwild .gt. 0)
     *               iseg(iwild(1)) = iseg(iwild(1)) + 1
                  do i = 2, nwild
                     if (iseg(iwild(i-1)) .le. 
     *                           eseg(iwild(i-1))) go to 100
                     iseg(iwild(i-1)) = bseg(iwild(i-1))
                     iseg(iwild(i))   = iseg(iwild(i)) + 1
                  enddo
  100             continue
               enddo
            else if (optable(c_opcode,iop) .eq. deallocate_op .or.
     *               (optable(c_opcode,iop) .eq. delete_op .and. 
     *                   partial_create) ) then

c---------------------------------------------------------------------------
c   Deallocate instruction: Wildcard indices are not encoded into the 
c   deallocate instruction, they must be decoded from the instruction that
c   last allocated the array.
c---------------------------------------------------------------------------

               array = optable(c_result_array,iop)
               instruction = allocate_table(array)
               allocate_table(array) = 0
               nblock = 1
               nwild  = 0
               nind   = 0
               do i = 1, mx_array_index
                  if (array_table(c_index_array1+i-1,array) .gt. 0) then
                     nind = nind + 1
                     ind(i) = array_table(c_index_array1+i-1,array)
                     nseg(i) = index_table(c_nsegments, ind(i))
                     bseg(i) = index_table(c_bseg, ind(i))
                     eseg(i) = index_table(c_eseg, ind(i))

                     if (optable(c_ind1+i-1,instruction) .eq. 
     *                                 wildcard_indicator) then
                        nblock = nblock*nseg(i)
                        nwild = nwild + 1
                        iwild(nwild) = i
                        iseg(iwild(nwild)) = bseg(iwild(nwild)) 
                     endif
                  endif
               enddo

               do j = 1, nind
                  maxrange(j) = 0
                  maxrange_seg(j) = 0
                  do jseg = bseg(j), eseg(j)
                     call get_index_segment(ind(j),jseg,
     *                             segment_table, nsegment_table,
     *                             index_table, nindex_table,
     *                             val1, val2)
                     if (val1 .le. 0 .and. val2 .le. 0) go to 900
                     n = val2 - val1 + 1  
                     if (n .gt. maxrange(j)) then
                        maxrange(j) = n
                        maxrange_seg(j) = jseg
                     endif
                  enddo
               enddo
               
               do iblock = 1, nblock
         
c--------------------------------------------------------------------------
c   Determine the blocksize.  First loop over the non-wildcard indices, 
c   using the maxrange for the size of each of their segments.
c--------------------------------------------------------------------------

                  n = 1
                  do i = 1, nind
                     if (optable(c_ind1+i-1,instruction) .ne. 
     *                            wildcard_indicator) then
                        n = n * maxrange(i)
                        seg(i) = maxrange_seg(i)
                     endif
                  enddo

c---------------------------------------------------------------------------
c   Next, loop over the wildcard indices, determining the actual segment
c   size for the next block.
c---------------------------------------------------------------------------

                  do i = 1, nwild
                     call get_index_segment(ind(iwild(i)),
     *                                      iseg(iwild(i)),
     *                             segment_table, nsegment_table,
     *                             index_table, nindex_table,
     *                             val1, val2)
                     n = n * (val2 - val1 + 1) 
                     seg(i) = iseg(iwild(i))
                  enddo

c--------------------------------------------------------------------------
c   Determine the stack on which to place the block.
c--------------------------------------------------------------------------

                  stack = which_stack(stack_blksizes, nstacks, n)

c---------------------------------------------------------------------------
c   Adjust the block count for the stack.
c---------------------------------------------------------------------------

                  if (partial_create) then

c---------------------------------------------------------------------------
c   Only delete it if the block_map_table indicates the block resides on
c   this processor.
c---------------------------------------------------------------------------

                     lookup =  block_map_lookup(seg, nind, array,
     *                       array_table(1,array),
     *                       index_table, nindex_table)
                     proc = block_map_table(c_processor,lookup)
                     if (proc .eq. my_company_rank) 
     *                   block_count(stack,iop) =
     *                        block_count(stack,iop) - 1
                  else
                     block_count(stack,iop) = 
     *                     block_count(stack,iop) - 1
                  endif

c--------------------------------------------------------------------------
c   Increment the "wildcard" segments.
c--------------------------------------------------------------------------

                  if (nwild .gt. 0)
     *               iseg(iwild(1)) = iseg(iwild(1)) + 1
                  do i = 2, nwild
                     if (iseg(iwild(i-1)) .le. 
     *                           eseg(iwild(i-1))) go to 200
                     iseg(iwild(i-1)) = bseg(iwild(i-1))
                     iseg(iwild(i))   = iseg(iwild(i)) + 1
                  enddo
  200             continue
               enddo
            endif

            do i = 1, nstacks
               running_count(i) = block_count(i,iop)
            enddo

  900    continue
         iop = iop + 1 
c            print *,'AT 900: iop, start_op, end_op = ',
c     *         iop, start_op, end_op
         if (start_op .ne. 0 .or. iop .le. noptable) then
            if (iop .gt. end_op)  iop = start_op
c            print *,'   BRANCH TO iop = ',iop
            if (iop .le. 0 .or. iop .gt. noptable) then
               print *,'ERROR: iop out of range: iop = ',iop
               call abort_job()
            endif
            go to 1000
         endif
         
 2000 continue

c---------------------------------------------------------------------------
c   Refine our block estimate by analyzing each loop, attempting to 
c   add blocks that are brought into and out of scope during loop execution.
c---------------------------------------------------------------------------

      do i = 1, max_level
         instr(i) = 0

         do j = 1, nstacks
            current(j,i) = 0
            running_count(j) = 0
         enddo
      enddo

      level = 1
      do iop = 1, noptable
         opcode = optable(c_opcode,iop)
         current_op = iop
         current_line = optable(c_lineno, iop)

         if (opcode .eq. pardo_op .or. 
     *       opcode .eq. do_op) then

c---------------------------------------------------------------------------
c   Increase the loop level and define the current number of blocks.
c---------------------------------------------------------------------------

            level = level + 1
            if (level .gt. max_level) then
               print *,'Error: Loop nesting is > ',max_level,
     *                 ' levels.'
               call abort_job()
            endif

            do j = 1, nstacks
               current(j,level) = 0
            enddo

            instr(level) = iop   ! save beginning instruction number

c--------------------------------------------------------------------------
c   The running count for the new loop is the sum of the newly
c   created blocks from all previous levels.
c--------------------------------------------------------------------------

            do j = 1, nstacks
               running_count(j) = 0
               do i = 1, level - 1
                  running_count(j) = 
     *                      running_count(j) + current(j,i)
               enddo
            enddo
         else if (opcode .eq. enddo_op .or. 
     *       opcode .eq. endpardo_op) then

c---------------------------------------------------------------------------
c   Move back one loop level.
c---------------------------------------------------------------------------

            level = level - 1
            if (level .lt. 1) then
               print *,'Error: Loop nesting level < 1 detected.'
               call abort_job()
            endif

c---------------------------------------------------------------------------
c   Compute a new running count.  
c---------------------------------------------------------------------------

            do j = 1, nstacks
               running_count(j) = 0
               do i = 1, level
                  running_count(j) = running_count(j) + current(j,i)
               enddo
            enddo
         else

c---------------------------------------------------------------------------
c   Bypass checking logic for logical instructions.
c---------------------------------------------------------------------------

            if (opcode .eq. jz_op .or.
     *          opcode .eq. go_to_op .or.
     *          opcode .eq. call_op .or.
     *          opcode .eq. destroy_op) go to 2300

            if (opcode .ge. sp_add_op .and.
     *          opcode .le. sp_ldindex_op) go to 2300

            if (opcode .ge. fl_add_op .and.
     *          opcode .le. fl_load_value_op) go to 2300

            if (opcode .eq. sp_ldi_sym_op) go to 2300

c---------------------------------------------------------------------------
c   Does this instruction create a temp result block?
c---------------------------------------------------------------------------

            array = optable(c_result_array,iop) 
            if (array .gt. narray_table) then
               print *,'ARRAY OUT OF BOUNDS: ',array,' iop ',
     *          iop,' op ',(optable(j,iop),j=1,loptable_entry)
               call abort_job()
            endif

            type  = array_table(c_array_type,array)

            if (opcode .eq. reindex_op) then
               do i = 1, mx_array_index
                  ind(i) = optable(c_ind1+i-1,iop)
               enddo

               call set_effective_indices(array_table(1,array),
     *                     ind)
            endif

            if (opcode .ne. reindex_op .and. array .ne. 0 ) then
               if (type .eq. temp_array .or.
     *             type .eq. served_array) then

c---------------------------------------------------------------------------
c   The current instruction results in creation of a result block.
c   Determine the maximum possible size of the block, and increment the
c   corresponding stack count.
c---------------------------------------------------------------------------

                  nind   = 0
                  do i = 1, mx_array_index
                     if (array_table(c_index_array1+i-1,array) .gt. 
     *                                             0) then
                        nind = nind + 1
                        ind(i) = array_table(c_index_array1+i-1,array)

                        if (ind(i) .gt. nindex_table) then
                           print *,'Task ',me,' INDEX OUT OF RANGE ',
     *                       ind(i),' iop ',iop,
     *                       ' op ',(optable(j,iop),j=1,loptable_entry)
                           print *,'Array ',array,' type ',type
                           call abort_job()
                        endif

                        nseg(i) = index_table(c_nsegments, ind(i))
                        bseg(i) = index_table(c_bseg, ind(i))
                        eseg(i) = index_table(c_eseg, ind(i))
                     endif
                  enddo

c---------------------------------------------------------------------------
c   Check all previous instructions at the current "level" to see if this
c   array and set of indices has already been used as a result.  If so, we 
c   need not count it.
c----------------------------------------------------------------------------

                  countit = .true.

                  do j = 1, nind
                     matchind(j) = 0
                  enddo

                  lopbegin = instr(level)
                  if (level .gt. 1) then
                     lopbegin = instr(2)  ! go through all nested loops.
                  endif

                  lopend = iop-1
                  if (optable(c_opcode,lopend) .eq. reindex_op .and.
     *                optable(c_result_array,lopend) .eq. array) 
     *               lopend = lopend - 1

                  do 2100 lop = lopbegin, lopend
                     if (optable(c_result_array,lop) .eq. array) then
                        if (optable(c_opcode,lop) .eq. reindex_op) then
                         
c--------------------------------------------------------------------------
c   Save indices from the "reindex" instruction.
c--------------------------------------------------------------------------

                           do j = 1, nind
                              matchind(j) = optable(c_ind1+j-1,lop)
                           enddo

                           go to 2100
                        endif

                        do j = 1, nind
                           if (matchind(j) .ne. ind(j)) then
                               go to 2100  ! try next instruction
                           endif
                        enddo 

                        countit = .false.   ! got a match, so do not count it
                        go to 2200
                     endif
 2100             continue
 2200             continue

                  if (countit) then

c----------------------------------------------------------------------------
c   Find the largest possible block with these indices.
c----------------------------------------------------------------------------

                     n = 1
                     do j = 1, nind
                        maxrange(j) = 0
                        do jseg = bseg(j), eseg(j)
                           call get_index_segment(ind(j),jseg,
     *                             segment_table, nsegment_table,
     *                             index_table, nindex_table,
     *                             val1, val2)
                           maxrange(j) = max(maxrange(j),(val2-val1+1))
                        enddo

                        n = n * maxrange(j)
                     enddo
               
c--------------------------------------------------------------------------
c   Determine the stack on which to place the block.
c--------------------------------------------------------------------------

                     stack = which_stack(stack_blksizes, nstacks, n)

c---------------------------------------------------------------------------
c   Adjust the block count for the stack.
c---------------------------------------------------------------------------

                     current(stack,level) = current(stack,level) + 1
                     running_count(stack) = running_count(stack) + 1      

c----------------------------------------------------------------------------
c   Add an additional block to account for prefetching the REQUEST instruction.
c------------------------------------------------------------------------------

                     if (opcode .eq. request_op) then
                        current(stack,level) = current(stack,level) + 1
                        running_count(stack) = running_count(stack) + 1
                     endif

                     block_count(stack,iop) =
     *                     block_count(stack,iop) + 
     *                        running_count(stack)
                  endif
               endif
            endif
 2300       continue
         endif
      enddo

c-------------------------------------------------------------------------
c   Restore array_table indices.
c-------------------------------------------------------------------------

      do i = 1, narray_table
         do j = 1, mx_array_index
            array_table(c_index_array1+j-1,i) = save_ind(j,i)
         enddo
      enddo

      call reset_timer_info()
      simulator = .false.
      return
      end

      integer function which_stack(stack_blksizes, nstacks, blocksize)
c------------------------------------------------------------------------
c   Determines the appropriate stack for a particular blocksize.
c   The stack_blksizes array is assumed to be ordered in increasing order.
c------------------------------------------------------------------------
      implicit none
      integer nstacks
      integer stack_blksizes(nstacks)
      integer blocksize
      integer i

      do i = 1, nstacks
         if (stack_blksizes(i) .ge. blocksize) then
            which_stack = i
            return
         endif
      enddo

      which_stack = nstacks   ! use largest blocksize.
      return
      end
