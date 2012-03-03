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
      subroutine pardo_loadb(optable, noptable, iop, index_table, 
     *                  nindex_table, array_table, narray_table,
     *                  block_map_table, segment_table, nsegment_table,
     *                  comm, debug, prt_flag,
     *                  start_op, end_op)
c---------------------------------------------------------------------------
c   PARDO loop with load balancing.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      INCLUDE 'timerz.h'
      include 'where_table.h'
      include 'batch_table.h'

      integer noptable, iop, start_op, end_op
      integer optable(loptable_entry,noptable)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer block_map_table(lblock_map_entry,*)
      integer nsegment_table 
      integer segment_table(lsegment_table_entry,nsegment_table) 
      integer comm

      integer i, j, opcode, istat
      integer result_array
      integer indx(mx_array_index)
      integer ind(mx_array_index), seg(mx_array_index)
      integer bseg(mx_array_index), eseg(mx_array_index)
      integer next_seg(mx_array_index)
      integer nind 
      integer my_lock
      integer pardo_master, get_pardo_master
      integer cluster_start, cluster_end
      logical debug, prt_flag
      logical prefetch

      integer push_do, pop_do

      integer next_batch, last_batch, max_batch

      save ind, nind, next_batch
      save bseg, eseg, seg

      if (iop .gt. noptable) then
	print *,'PARDO_LOOP: Called with iop = ',iop,' noptable = ',
     *           noptable
        call abort_job()
      endif

      if (iop .lt. start_op .or. iop .gt. end_op) return
      opcode = optable(c_opcode, iop)

      if (opcode .eq. pardo_op) then
         pardo_master = get_pardo_master()

c--------------------------------------------------------------------------
c   PARDO operation
c   Save the loop indices, line number.
c--------------------------------------------------------------------------

         do i = 1, mx_array_index
            indx(i) = optable(c_ind1+i-1,iop)
         enddo

         call set_prefetch_context(indx, mx_array_index)
 
c---------------------------------------------------------------------------
c   First pass through the loop?
c---------------------------------------------------------------------------

         my_lock = optable(c_pardo_lock_index,iop)
         if (my_company_rank .eq. pardo_master) 
     *          call f_acquire_pardo_lock(my_lock)
         if (optable(c_oploop, iop) .eq. 0) then
c            print *,'TASK ',me,' PARDO INIT CODE line ',current_line
c            call c_flush_stdout()

c----------------------------------------------------------------------------
c   Push the start_op and end_op indices on the instruction stack.
c----------------------------------------------------------------------------

            istat = push_do(start_op, end_op)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *           ' Error in pardo_loadb: stack overflow'
               call abort_job()
            endif

c---------------------------------------------------------------------------
c   Push the data for the WHERE conditionals.
c---------------------------------------------------------------------------

            istat = push_do(iwhere, nwhere)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *           ' Error in pardo_loadb: stack overflow while ',
     *           'pushing wheres'
               call abort_job()
            endif

c---------------------------------------------------------------------------
c   Set the new values of start_op, end_op
c---------------------------------------------------------------------------

            start_op = iop
            do j = iop+1,noptable
               if (optable(c_opcode,j) .eq. endpardo_op) then
                  end_op = j
                  go to 100
               endif
            enddo
  100       continue

            optable(c_oploop, iop)     = 1   ! loop is initialized

c--------------------------------------------------------------------------
c   Get the current pardo timer keys out of the optable.
c--------------------------------------------------------------------------

            call unpack_pardo_timer(optable(c_instr_timer,iop),
     *                      pardo_timer, pardo_block_wait_timer) 
            call timer_start(pardo_timer)

c----------------------------------------------------------------------------
c   Get the WHERE conditionals for the current loop.
c----------------------------------------------------------------------------

            call get_where_conditionals(optable, start_op, end_op)

c----------------------------------------------------------------------------
c   Set up loop segment/index mapping.
c----------------------------------------------------------------------------

            max_batch = -1  
            call  lb_set_loop_incr_mapping(indx, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      ind, bseg, eseg, nind,
     *                      max_batch)
            last_batch_processed = optable(c_pardo_max_batch,iop) + 1
         endif   

c--------------------------------------------------------------------------
c   Determine the next loop batch.
c--------------------------------------------------------------------------

         next_batch = optable(c_pardo_batch,iop) + 1
         last_batch = optable(c_pardo_batch_end,iop)

         if (optable(c_pardo_batch,iop) .ge. 0) 
     *      optable(c_pardo_batch,iop) = optable(c_pardo_batch,iop)+1  

         if (next_batch .gt. last_batch .or.
     *       next_batch .eq. 0) then
            if (my_company_rank .eq. pardo_master) then
               call pardo_loadb_update_batch(iop, optable, index_table, 
     *                segment_table, nsegment_table,
     *                next_batch, last_batch)

               call exec_thread_server(0) 
               if (next_batch .le. 0) optable(c_oploop,iop) = 0
            else
               call f_acquire_pardo_lock(my_lock)
               call pardo_loadb_get_next_batch(comm, iop, index_table,
     *                 next_batch, last_batch)

               optable(c_pardo_batch, iop) = next_batch
               optable(c_pardo_batch_end,iop) = last_batch

               if (next_batch .le. 0) optable(c_oploop,iop) = 0
               call f_release_pardo_lock(my_lock)
            endif   
         endif

         if (my_company_rank .eq. pardo_master) 
     *           call f_release_pardo_lock(my_lock)

c---------------------------------------------------------------------------
c   We now have a batch.  Check for end of loop.
c---------------------------------------------------------------------------

         if (next_batch .le. 0) then
            last_batch_processed = next_batch
            iop = end_op                ! point instruction counter to endpardo
            return                      ! Next instruction will be the endpardo
         endif

         call lb_get_next_loop_set(next_batch, ind, nind,
     *                             index_table, nindex_table,
     *                             bseg, eseg, seg)

c--------------------------------------------------------------------------
c   Save batch info in the batch table.
c--------------------------------------------------------------------------
 
         last_batch_processed = next_batch
         do i = 1, nind
            save_batch(i) = seg(i)
         enddo

         call lb_get_next_loop_set(next_batch+1, ind, nind,
     *                             index_table, nindex_table,
     *                             bseg, eseg, next_seg)

c---------------------------------------------------------------------------
c   Store the set of index segments in the "current_seg" field of the 
c   index table.
c---------------------------------------------------------------------------

         prefetch = .true.
         do i = 1, nind
            index_table(c_current_seg,ind(i)) = seg(i)
            index_table(c_next_seg, ind(i))   = next_seg(i)
            if (next_seg(i) .eq. 0) prefetch = .false.
         enddo

c----------------------------------------------------------------------------
c   If we cannot prefetch, zero out all the "next_seg" entries.
c----------------------------------------------------------------------------

         if (.not. prefetch) then
            do i = 1, nind
               index_table(c_next_seg, ind(i)) = 0
            enddo
         endif 
         
         iop = start_op + 1   ! point to next instruction
         return
      endif   ! pardo_op

      if (opcode .eq. endpardo_op) then
c         print *,'Task ',me,' Endpardo at line ',current_line
c         call c_flush_stdout()
         call unset_prefetch_context()

c-------------------------------------------------------------------------
c   Reset any temp blocks created in this loop.
c-------------------------------------------------------------------------

         do i = start_op, end_op

c-------------------------------------------------------------------------
c   If a PUT was executed, set the flag in the array_table.
c-------------------------------------------------------------------------
                                                                                
            if (optable(c_opcode,i) .eq. put_op .or.
     *          optable(c_opcode,i) .eq. put_replace_op) then
               result_array = optable(c_result_array,i)
               array_table(c_put_flag,result_array) = 1
            endif

c-------------------------------------------------------------------------
c   If a PREPARE was executed, set the flag in the array_table.
c-------------------------------------------------------------------------
                                                                                
            if (optable(c_opcode,i) .eq. prepare_op .or.
     *          optable(c_opcode,i) .eq. prepare_increment_op) then
               result_array = optable(c_result_array,i)
               array_table(c_prepare_flag,result_array) = 1
            endif

            if (optable(c_opblock,i) .ne. 0) then

c-------------------------------------------------------------------------
c   Only clear the block_computed_flag if the array is managed by
c   blkmgr.
c-------------------------------------------------------------------------

               result_array = optable(c_result_array,i)
               if (array_table(c_array_type, result_array) .ne.
     *             static_array) then
                  call block_end_of_loop(result_array,
     *              optable(c_opblock,i),optable(c_opblkndx,i),
     *              array_table, narray_table, index_table,
     *              nindex_table, block_map_table)
               endif
               optable(c_opblock,i) = 0
               optable(c_opblkndx,i) = 0
            endif
         enddo

c-------------------------------------------------------------------------
c   ENDPARDO operation.
c   Check for loop termination.
c-------------------------------------------------------------------------

         if (next_batch .lt. 0) then 
            optable(c_pardo_batch,iop) = -1
            call update_timer(pardo_timer)

c            print *,'Task ',me,' ENDPARDO loop termination'

c--------------------------------------------------------------------------
c   Reset the loop indices to "undefined".
c--------------------------------------------------------------------------

            do i = 1, mx_array_index
               ind(i) = optable(c_ind1+i-1,iop)
               if (ind(i) .ge. 1 .and. ind(i) .le. nindex_table) then
                  index_table(c_current_seg,ind(i)) = undefined_segment
                  index_table(c_next_seg, ind(i))   = 0
               endif  
            enddo

c--------------------------------------------------------------------------
c   Final pass of the loop.  Pop a new set of iwhere, nwhere values off the
c   stack.
c--------------------------------------------------------------------------

            istat = pop_do(iwhere, nwhere)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *          ' Error in pardo_loadb: stack underflow'
               call abort_job()
            endif

c--------------------------------------------------------------------------
c   Final pass of the loop.  Pop a new set of start, end values off the
c   stack.
c--------------------------------------------------------------------------

            istat = pop_do(start_op, end_op)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *          ' Error in pardo_loadb: stack underflow'
               call abort_job()
            endif

            if (start_op .gt. noptable .or. 
     *          start_op .lt. 0 .or.
     *          end_op .gt. noptable .or.
     *          end_op .lt. 0) call abort_job()
         endif 

      endif   ! endpardo_op

      return
      end

      subroutine set_pardo_workload(max_batch, start_batch, end_batch)
c--------------------------------------------------------------------------
c   Defines the pardo workload for the current pardo cluster. 
c--------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
      include 'trace.h'
      integer max_batch, start_batch, end_batch
      integer get_pardo_master, pardo_master
      integer cluster_size, get_pardo_cluster_size
      integer nbatch, ndiff, nclusters, my_cluster

      cluster_size = get_pardo_cluster_size()
      pardo_master = get_pardo_master()
      nclusters     = my_company_size / cluster_size
      if (nclusters * cluster_size .lt. my_company_size)
     *   nclusters = nclusters + 1
      nbatch       = max_batch / nclusters
      ndiff        = max_batch - nbatch * nclusters
      my_cluster   = pardo_master / cluster_size 
      start_batch  = my_cluster * nbatch + 1
      if (ndiff .gt. 0) then
         if (my_cluster .lt. ndiff) then
            nbatch = nbatch + 1
            start_batch = start_batch + my_cluster  
         else
            start_batch = start_batch + ndiff
         endif
      endif

      end_batch = start_batch + nbatch - 1 
      return
      end

      subroutine lb_set_loop_incr_mapping(raw_indx, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      ind, bseg, eseg, nind,
     *                      max_batch)
c--------------------------------------------------------------------------
c   Defines a new parallel mapping of the loop indices in array "ind", 
c   based on the index_table.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'where_table.h'
      include 'batch_table.h'
      include 'parallel_info.h'
      include 'trace.h'

      integer raw_indx(mx_array_index)
      integer indx(mx_array_index)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer nsegment_table
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer i, j, nind, bseg(mx_array_index),eseg(mx_array_index)
      integer ind(mx_array_index)
      integer max_batch
      integer superseg, superindex

      integer seg(mx_array_index)
      integer count, ncomb
      logical check_where_conditions
      logical subrange_flag(mx_array_index)

c--------------------------------------------------------------------------
c   Count the number of indices, determine number of segments.
c--------------------------------------------------------------------------

      nind = 0
      do i = 1, mx_array_index
         if (raw_indx(i) .ne. 0) then
            call get_loop_index(raw_indx(i), indx(i), subrange_flag(i))

            nind = nind + 1
            ind(nind) = indx(i)

            if (subrange_flag(i)) then

c---------------------------------------------------------------------------
c   Find the superindex and its current segment (superseg).
c---------------------------------------------------------------------------

               superindex = index_table(c_subindex_ptr,ind(nind))
               superseg   = index_table(c_current_seg,superindex)
               if (superseg .eq. undefined_segment) then
                  print *,'Error: Pardo at line ',current_line,
     *              ' superindex ',superindex,' is not in scope.' 
                  call abort_job()
               endif
  
c--------------------------------------------------------------------------
c   Find the segment_table entry for superseg.
c--------------------------------------------------------------------------

               do j = 1, nsegment_table
                  if (segment_table(c_index,j) .eq. superindex .and.
     *                segment_table(c_segment,j) .eq. superseg) then
                     bseg(nind) = segment_table(c_subseg1,j)
                     eseg(nind) = segment_table(c_subseg2,j)
                     go to 50
                  endif 
               enddo
   50          continue
            else
               eseg(nind) = index_table(c_eseg, indx(i))
               bseg(nind) = index_table(c_bseg, indx(i))
            endif
         endif
      enddo

      if (max_batch .eq. -1) return

c--------------------------------------------------------------------------
c   Determine the max. number of batches to process.
c--------------------------------------------------------------------------

      ncomb = 1
      do i = 1, nind
         ncomb = ncomb * (eseg(i)-bseg(i)+1)
      enddo

      if (nwhere .eq. 0) then

c---------------------------------------------------------------------------
c   No WHERE conditions are applied, simply multiply the ranges of all the 
c   indices together.
c---------------------------------------------------------------------------

         max_batch = ncomb
      else

c---------------------------------------------------------------------------
c   When we have WHERE conditions, we must count the number of valid 
c   combinations of indices.
c---------------------------------------------------------------------------

         do i = 2, nind
            seg(i) = bseg(i)
         enddo
         seg(1) = bseg(1) - 1

         max_batch = 0
         do count = 1, ncomb
            do i = 1, nind
               seg(i) = seg(i) + 1
               if (seg(i) .le. eseg(i)) then
                  go to 100
               else
                  seg(i) = bseg(i)
               endif
            enddo

  100       continue

c---------------------------------------------------------------------------
c   Check the combination for valid WHERE conditions.
c---------------------------------------------------------------------------

            if (check_where_conditions(ind, seg, nind, 
     *                                 index_table, nindex_table)) 
     *            max_batch = max_batch + 1
         enddo
      endif

c---------------------------------------------------------------------------
c   Initialize the batch_table.
c---------------------------------------------------------------------------

      last_batch_processed = max_batch + 1   ! force start at 1.
      return
      end

      subroutine lb_get_next_loop_set(next_batch, ind, nind,
     *                                index_table, nindex_table,
     *                                bseg, eseg, seg)
      implicit none
      include 'maxdim.h'
      include 'where_table.h'
      include 'batch_table.h'
 
      integer next_batch, nind
      integer nindex_table
      integer index_table(*)
      integer ind(nind), bseg(nind), eseg(nind)
      integer seg(nind)

      integer i
      integer count
      integer istart
      
      logical check_where_conditions

c--------------------------------------------------------------------------
c   Convert the batch index into segment indices.
c--------------------------------------------------------------------------

      if (last_batch_processed .gt. next_batch .or.
     *    last_batch_processed .le. 0) then

c--------------------------------------------------------------------------
c   Start at batch number 1.
c--------------------------------------------------------------------------

         istart = 1
         do i = 2, nind
            seg(i) = bseg(i)
         enddo
         seg(1) = bseg(1) - 1
      else

c---------------------------------------------------------------------------
c   Start at batch number 'last_batch'
c---------------------------------------------------------------------------

         istart = last_batch_processed+1
         do i = 1, nind
            seg(i) = save_batch(i)
         enddo  
      endif

      do count = istart, next_batch
   50    continue      
         do i = 1, nind
            seg(i) = seg(i) + 1 
            if (seg(i) .le. eseg(i)) then
               go to 100
            else
               seg(i) = bseg(i) 
            endif
         enddo

  100    continue

c--------------------------------------------------------------------------
c   If there are any "where" conditions, make sure they are satisfied.
c--------------------------------------------------------------------------

         if (nwhere .gt. 0) then
            if (.not. check_where_conditions(ind, seg, nind,
     *                                 index_table, nindex_table))
     *            go to 50
         endif

      enddo 

      return
      end
