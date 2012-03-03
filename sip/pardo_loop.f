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
      subroutine pardo_loop(optable, noptable, iop, index_table, 
     *                  nindex_table, array_table, narray_table,
     *                  block_map_table,
     *                  segment_table, nsegment_table,
     *                  comm, debug, prt_flag,
     *                  start_op, end_op)

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'where_table.h'

      integer iblk

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

      integer i, j, opcode, index, nseg, istat
      integer result_array
      integer ind(mx_array_index), segm(mx_array_index)
      integer indx(mx_array_index), bseg(mx_array_index),
     *        eseg(mx_array_index), max_batch, nind
      integer next_seg(mx_array_index)
      integer first_batch, last_batch, current_batch
      integer nsegm
      logical debug, prt_flag

      integer push_do, pop_do

      save indx, bseg, eseg, max_batch, nind 
      save first_batch, last_batch, current_batch

      opcode = optable(c_opcode, iop)
      if (iop .gt. noptable) then
         print *,'PARDO_LOOP: Called with iop = ',iop,' noptable = ',
     *           noptable
         call abort_job()
      endif

      if (iop .lt. start_op .or. iop .gt. end_op) return
      if (opcode .eq. pardo_op) then
c         print *,'PARDO at line ',current_line

c--------------------------------------------------------------------------
c   PARDO operation
c   Update the "current_seg" fields of all loop indices.
c--------------------------------------------------------------------------

         do i = 1, mx_array_index
            ind(i) = optable(c_ind1+i-1,iop)
         enddo

         call set_prefetch_context(ind, mx_array_index)

         if (optable(c_oploop,iop) .eq. 1) then

c---------------------------------------------------------------------------
c   Find the segments corresponding to the next "batch" of work to process.
c---------------------------------------------------------------------------
 
            current_batch = current_batch + 1
            call get_current_loop_set(current_batch, first_batch,
     *                                last_batch, indx, bseg, eseg,
     *                                nind, index_table, nindex_table, 
     *                                segm)

c---------------------------------------------------------------------------
c   Store the set of index segments in the "current_seg" field of the 
c   index table.
c---------------------------------------------------------------------------

            do i = 1, nind
               index = indx(i)
               if (index .gt. 0) 
     *             index_table(c_current_seg,index) = segm(i)
            enddo

c---------------------------------------------------------------------------
c   Check to see if we have any work to do?
c---------------------------------------------------------------------------

            nsegm = 0
            do i = 1,nind
               if (segm(i) .gt. 0) nsegm = nsegm + 1
            enddo
c            print *,'Pardo new batch segm ',(segm(i),i=1,mx_array_index)

c            print *,'Task ',me,' current_batch ',current_batch,
c     *        ' segm ',(segm(i),i=1,nind)

            if (nsegm .gt. 0) then

c---------------------------------------------------------------------------
c   We have a legitimate set of segments for processing.  Now look at the
c   next batch and set up the next_seg field of the index_table as 
c   prefetch segments.
c---------------------------------------------------------------------------

               call get_current_loop_set(current_batch+1, first_batch,
     *                                last_batch, indx, bseg, eseg,
     *                                nind, index_table, nindex_table, 
     *                                next_seg)
               do i = 1, nind
                  index_table(c_next_seg, indx(i)) = next_seg(i)
               enddo
            else
               print *,'Error: PARDO at line ',current_line,
     *             ' has no valid segments to process.'
               print *,'   Task me = ',me,' segm = ',
     *                     (segm(i),i=1,mx_array_index)
               call abort_job()
            endif
         else

c----------------------------------------------------------------------------
c   LOOP INITIALIZATION:
c   Push the start_op and end_op indices on the instruction stack.
c----------------------------------------------------------------------------

            istat = push_do(start_op, end_op)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *           ' Error in pardo_loop: stack overflow'
               call abort_job()
            endif

c---------------------------------------------------------------------------
c   Push the data for the WHERE conditionals.
c---------------------------------------------------------------------------

            istat = push_do(iwhere, nwhere)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *           ' Error in pardo_loop: stack overflow while ',
     *           'pushing wheres'
               call abort_job()
            endif

c---------------------------------------------------------------------------
c   Unpack the current pardo timers and start the timer.
c---------------------------------------------------------------------------

            call unpack_pardo_timer(optable(c_instr_timer,iop),
     *                 pardo_timer, pardo_block_wait_timer)
            call timer_start(pardo_timer)

c---------------------------------------------------------------------------
c   Set the new values of start_op, end_op
c---------------------------------------------------------------------------

c            print *,'Old start_op, end_op ',start_op,end_op 
            start_op = iop
            do j = iop+1,noptable
               if (optable(c_opcode,j) .eq. endpardo_op) then
                  end_op = j
                  go to 100
               endif
            enddo
  100       continue
c            print *,'Task ',me,' PARDO_INIT code: new iop = ',
c     *              iop,' start_op, end_op = ',start_op,end_op

c----------------------------------------------------------------------------
c   Get the WHERE conditionals for the current loop.
c----------------------------------------------------------------------------

            call get_where_conditionals(optable, start_op, end_op)

c---------------------------------------------------------------------------
c   Entering loop for the first iteration: set up the mapping for 
c   incrementing the indices, then pull out the first index set.
c---------------------------------------------------------------------------

            optable(c_oploop,iop) = 1

            max_batch = -1 
            call lb_set_loop_incr_mapping(ind, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      indx, bseg, eseg, nind,
     *                      max_batch)
            max_batch = optable(c_pardo_max_batch,iop)
            call set_loop_incr_mapping(first_batch, last_batch,
     *                                 max_batch) 
            current_batch = first_batch
            call get_current_loop_set(current_batch, first_batch, 
     *                                last_batch, indx, bseg, eseg, 
     *                                nind, index_table, nindex_table, 
     *                                segm)

c---------------------------------------------------------------------------
c   Store the set of index segments in the "current_seg" field of the 
c   index table.
c---------------------------------------------------------------------------

            nsegm = 0
            do i = 1, nind
               if (segm(i) .gt. 0) then
                  nsegm = nsegm + 1
                  index = indx(i)
                  if (index .gt. 0) 
     *                index_table(c_current_seg,index) = segm(i)
               endif
            enddo

c---------------------------------------------------------------------------
c   If we got a real set of segments to process, set the "next_seg"
c   field used for prefetching data.
c---------------------------------------------------------------------------

            if (nsegm .gt. 0) then
               call get_current_loop_set(current_batch+1, first_batch,
     *                                last_batch, indx, bseg, eseg,
     *                                nind, index_table, nindex_table,
     *                                next_seg)
               do i = 1, nind
                  index_table(c_next_seg, indx(i)) = next_seg(i)
               enddo
            else
c               if (me .eq. 0) print *,'PARDO line ',current_line,
c     *              ' END OF BATCH'
               do i = 1, nind
                  if (indx(i) .gt. 0) 
     *               index_table(c_next_seg, indx(i)) = 0
               enddo

c---------------------------------------------------------------------------
c   No work for this processor.  Set iop to point to the "endpardo" 
c   instruction and return.  This is equivalent to an unconditional 
c   jump to the end of the loop.
c---------------------------------------------------------------------------

               optable(c_oploop,iop) = 0   ! reset loop init flag
               current_batch = last_batch
               iop = end_op
               return
            endif   ! nsegm .gt. 0
         endif      ! pardo loop init code
      endif         ! pardo_op

      if (opcode .eq. endpardo_op) then
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

c-------------------------------------------------------------------------
c   Only reset the block's flags if this is NOT the result of a CREATE
c   instruction.  Only a DELETE can reset the flags otherwise.
c-------------------------------------------------------------------------

                 call block_end_of_loop(result_array,
     *                 optable(c_opblock,i), optable(c_opblkndx,i),
     *                 array_table, narray_table, index_table,
     *                 nindex_table, block_map_table)
               endif
               optable(c_opblock,i) = 0
               optable(c_opblkndx,i) = 0
            endif
         enddo

         if (current_batch .eq. last_batch) then
c
c-----------------------------------------------------------------------
c   Reset the loop indices to "undefined".
c-----------------------------------------------------------------------

         do i = 1, nind
            if (indx(i) .ge. 1 .and. 
     *          indx(i) .le. nindex_table) then
               index_table(c_current_seg,indx(i)) = undefined_segment
               index_table(c_next_seg,indx(i))    = 0
            endif
         enddo

c--------------------------------------------------------------------------
c   Final pass of the loop.  Pop a new set of iwhere, nwhere values off the
c   stack.
c--------------------------------------------------------------------------

            istat = pop_do(iwhere, nwhere)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *          ' Error in pardo_loop: stack underflow'
               call abort_job()
            endif

c--------------------------------------------------------------------------
c   Final pass of the loop.  Turn off the pardo_timer and pop a new set 
c   of start, end values off the stack.
c--------------------------------------------------------------------------

            optable(c_oploop,start_op) = 0   ! reset loop init flag
            call update_timer(pardo_timer)

            istat = pop_do(start_op, end_op)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *          ' Error in pardo_loop: stack underflow'
               call abort_job()
            endif
     
            if (start_op .gt. noptable .or. 
     *          start_op .lt. 0 .or.
     *          end_op .gt. noptable .or.
     *          end_op .lt. 0) call abort_job()
         endif
      endif

      return
      end

      subroutine set_loop_incr_mapping(first_batch, last_batch, 
     *                                 nbatch_total)
c--------------------------------------------------------------------------
c   Defines a new parallel mapping of the loop indices in array "ind", 
c   based on the index_table.
c--------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'parallel_info.h'
      include 'trace.h'
      include 'where_table.h'

      integer ierr

      integer i
      integer my_proc, first_batch, last_batch, nbatch
      integer nbatch_total, nbatch_proc, nbatch_left
      integer next_batch, current_batch

      my_proc = my_company_rank

      nbatch_proc = nbatch_total / my_company_size
      nbatch_left = nbatch_total - my_company_size * nbatch_proc

      if (my_proc .lt. nbatch_left) then

c-------------------------------------------------------------------------
c   The first "nbatch_left" procs have 1 extra batch to do.
c-------------------------------------------------------------------------

         first_batch = my_proc * (nbatch_proc + 1) + 1
         last_batch  = first_batch + nbatch_proc
      else

c--------------------------------------------------------------------------
c   Procs "nbatch_left" to "my_company_size" have only nbatch_proc batches 
c   to do.
c--------------------------------------------------------------------------

         first_batch = my_proc * nbatch_proc + 
     *                  nbatch_left + 1
         last_batch = first_batch + nbatch_proc - 1
      endif
          
      nbatch = last_batch - first_batch + 1
      return
      end

      subroutine get_current_loop_set(current_batch, first_batch,
     *                                last_batch, ind, bseg, eseg,
     *                                nind, index_table, nindex_table,
     *                                segm)
      implicit none
      include 'where_table.h'
      
      integer current_batch, first_batch, last_batch
      integer i, ibatch, nind
      integer ind(nind), bseg(nind), eseg(nind), segm(nind)
      integer index_table, nindex_table

      logical check_where_conditions
 
c--------------------------------------------------------------------------
c   Breaks down the next set of loop coordinates, returns them in "seg".
c   If the current_batch > last_batch, the segm argument is returned
c   with all 0's.  Otherwise the segments of the batch are returned in 
c   "segm".
c--------------------------------------------------------------------------

      if (current_batch > last_batch) then
         do i = 1, nind
            segm(i) = 0
         enddo
 
         return
      endif

c---------------------------------------------------------------------------
c  Initialize the segm array to the bseg for each index.
c---------------------------------------------------------------------------

      do i = 1, nind
         segm(i) = bseg(i)
      enddo
      segm(1) = segm(1) - 1

c---------------------------------------------------------------------------
c   Cycle the indices to "current_batch".
c---------------------------------------------------------------------------

      do ibatch = 1, current_batch
   50    continue
         do i = 1, nind
            segm(i) = segm(i) + 1
            if (segm(i) .le. eseg(i)) then
               go to 100
            else
               segm(i) = bseg(i)
            endif
         enddo
  100    continue

c--------------------------------------------------------------------------
c   If there are any "where" conditions, make sure they are satisfied.
c--------------------------------------------------------------------------

         if (nwhere .gt. 0) then
            if (.not. check_where_conditions(ind, segm, nind,
     *                                 index_table, nindex_table))
     *            go to 50
         endif
      enddo

      return
      end
