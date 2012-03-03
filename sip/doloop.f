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
      subroutine doloop(optable, noptable, iop, index_table, 
     *                  nindex_table, array_table, narray_table,
     *                  block_map_table, segment_table, nsegment_table,
     *                  debug, prt_flag,
     *                  start_op, end_op)

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      include 'where_table.h'

      integer noptable, iop, start_op, end_op
      integer optable(loptable_entry,noptable)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer block_map_table(lblock_map_entry,*)
      integer nsegment_table
      integer segment_table(lsegment_table_entry,nsegment_table)

      integer i, j, opcode, index, iseg, istat
      integer result_array
      integer bseg, final_seg, next_seg
      integer superindex, superseg 
      logical debug, prt_flag
      logical check_where_conditions
      logical subrange_flag

      integer push_do, pop_do


      opcode = optable(c_opcode, iop)
      if (iop .gt. noptable) then
         print *,'DOLOOP: Called with iop = ',iop,' noptable = ',
     *           noptable
         call abort_job()
      endif

      if (iop .lt. start_op .or. iop .gt. end_op) return
      if (opcode .eq. do_op) then

c--------------------------------------------------------------------------
c   DO operation
c   Increment the "current_seg" field of the loop index.
c--------------------------------------------------------------------------

         call get_loop_index(optable(c_ind1,iop), index, subrange_flag)
         call set_prefetch_context(index, 1)
         if (optable(c_oploop,iop) .eq. 1) then
            index_table(c_current_seg,index) = 
     *        index_table(c_current_seg,index) + 1

c---------------------------------------------------------------------------
c   Update the c_next_seg entry (used in prefetches).
c---------------------------------------------------------------------------

            next_seg = index_table(c_current_seg,index) + 1
            if (next_seg .gt. index_table(c_eseg,index)) next_seg = 0 
            index_table(c_next_seg,index) = next_seg

            if (prt_flag) 
     *         print *,'Task ',me,' index = ',index,' current_seg = ',
     *                 index_table(c_current_seg,index)
            if (debug)
     *         print *,'DO index ',index,' set to ',
     *            index_table(c_current_seg,index)

         else

c--------------------------------------------------------------------------
c   New loop initialization.
c--------------------------------------------------------------------------
         
            call doloop_init(optable, noptable, iop, index_table,
     *                  nindex_table, segment_table, nsegment_table,
     *                  debug, prt_flag,
     *                  start_op, end_op, index, subrange_flag)
         endif

         if (nwhere .gt. 0) then

c---------------------------------------------------------------------------
c   If the segment fails the WHERE conditions, go to end of loop.
c---------------------------------------------------------------------------

            if (.not. check_where_conditions(index, 
     *                       index_table(c_current_seg,index), 1, 
     *                       index_table, nindex_table) ) then 
               iop = end_op
               return 
            endif
         endif
      endif

      if (opcode .eq. enddo_op) then
         call unset_prefetch_context()   ! back out the current context.

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
c   Only clear the block_computed_flag if the block is managed by
c   blkmgr.
c-------------------------------------------------------------------------

               result_array = optable(c_result_array,i)
               if (array_table(c_array_type, result_array) .ne. 
     *             static_array) then
                  call block_end_of_loop(result_array, 
     *                  optable(c_opblock,i), optable(c_opblkndx,i),
     *                  array_table, narray_table, index_table,
     *                  nindex_table, block_map_table)
               endif   
               optable(c_opblock,i) = 0
               optable(c_opblkndx,i) = 0
            endif   
         enddo

c-------------------------------------------------------------------------
c   ENDDO operation.
c   Check for final iteration of loop.
c-------------------------------------------------------------------------

         call get_loop_index(optable(c_ind1,iop), index, subrange_flag)
         if (index .ge. 1 .and. index .le. nindex_table) then
            iseg  = index_table(c_current_seg,index)
            if (debug)
     *      print *,'ENDDO for index ',index,' iseg, current_seg = ',
     *         iseg,index_table(c_nsegments,index) 

            if (subrange_flag) then

c--------------------------------------------------------------------------
c   Look up the superindex of this subindex.
c--------------------------------------------------------------------------

               superindex = index_table(c_subindex_ptr,index)
               superseg   = index_table(c_current_seg,superindex)
               if (superseg .eq. undefined_segment) then
                  print *,'Error: Superindex ',superindex,
     *              ' not in scope at line ',current_line
                  call abort_job()
               endif 
   
c--------------------------------------------------------------------------
c   Find the c_subseg2 field of the current segment of the superindex.
c--------------------------------------------------------------------------

               do i = 1, nsegment_table
                  if (segment_table(c_index,i) .eq. superindex .and.
     *                segment_table(c_segment,i) .eq. superseg) then
                     final_seg = segment_table(c_subseg2,i)
                     go to 100
                  endif
               enddo   
  100          continue  
            else
               final_seg = index_table(c_eseg,index)
            endif

            if (iseg .eq. final_seg) then

c----------------------------------------------------------------------------
c   Reset the loop init flag.
c----------------------------------------------------------------------------

               optable(c_oploop, start_op) = 0   ! should be a do instruction

c--------------------------------------------------------------------------
c   Reset the loop index segment to "undefined".
c--------------------------------------------------------------------------

               index_table(c_current_seg, index) = undefined_segment
               index_table(c_next_seg, index) = 0

c--------------------------------------------------------------------------
c   Final pass of the loop.  Pop a new set of iwhere, nwhere values off the
c   stack.
c--------------------------------------------------------------------------

               istat = pop_do(iwhere, nwhere)
               if (istat .lt. 0) then
                  print *,'Task ',me,
     *             ' Error in doloop: stack underflow'
                  call abort_job()
               endif

c--------------------------------------------------------------------------
c   Final pass of the loop.  Pop a new set of start, end values off the
c   stack.
c--------------------------------------------------------------------------

               istat = pop_do(start_op, end_op)
               if (istat .lt. 0) then
                  print *,'Task ',me,
     *             ' Error in doloop: stack underflow'
                  call abort_job()
               endif
            endif
         endif
      endif

      return
      end
       
      subroutine doloop_init(optable, noptable, iop, index_table,
     *                  nindex_table, segment_table, nsegment_table,
     *                  debug, prt_flag,
     *                  start_op, end_op, loop_index, subrange_flag)
c----------------------------------------------------------------------------
c   Initialization code for do loop.
c   Called upon entry to a new loop.
c----------------------------------------------------------------------------
      implicit none

      include 'interpreter.h'
      include 'where_table.h'
      include 'parallel_info.h'
      include 'trace.h'

      integer noptable, iop, start_op, end_op
      integer optable(loptable_entry,noptable)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer nsegment_table
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer loop_index, instr_index
      integer bseg, eseg, next_seg
      integer i, superindex, superseg
      logical debug, prt_flag
      logical subrange_flag
      logical dummy

      integer opcode, index, istat, j
      integer push_do

      if (optable(c_oploop, iop) .ne. 0) return   ! loop init is already done

      opcode = optable(c_opcode, iop)
      if (opcode .ne. do_op) return

      index = loop_index
      if (subrange_flag) then

c---------------------------------------------------------------------------
c   We are looping over the blocks of a subindex that are contained in a
c   single superindex segment.  First look up the superindex.
c---------------------------------------------------------------------------

         superindex = index_table(c_subindex_ptr,index)
        
c---------------------------------------------------------------------------
c   Make sure the superindex is in scope.
c---------------------------------------------------------------------------
   
         superseg = index_table(c_current_seg,superindex)
         if (superseg .eq. undefined_segment) then
            print *,'Error: Superindex ',superindex,' not in scope at ',
     *                'line ',current_line
            call abort_job() 
         endif

c---------------------------------------------------------------------------
c   Set the current segment of the subindex to the c_subseg1 entry of its
c   superindex.
c---------------------------------------------------------------------------
   
         do i = 1, nsegment_table
            if (segment_table(c_index,i) .eq. superindex .and.
     *          segment_table(c_segment,i) .eq. superseg) then
               index_table(c_current_seg,index) =
     *                     segment_table(c_subseg1,i)
               go to 100
            endif 
         enddo
  100    continue
      else
         index_table(c_current_seg,index) = index_table(c_bseg, index)
         next_seg = index_table(c_current_seg,index) + 1
         if (next_seg .gt. index_table(c_eseg,index)) next_seg = 0
         index_table(c_next_seg,index) = next_seg 
      endif

      if (debug) 
     *   print *,'NEW DO index ',index,' set to segment ',
     *     index_table(c_current_seg,index)

      istat = push_do(start_op, end_op)
      if (istat .lt. 0) then
         print *,' Error in doloop: stack overflow'
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Push the data for the WHERE conditionals.
c---------------------------------------------------------------------------

      istat = push_do(iwhere, nwhere)
      if (istat .lt. 0) then
         print *,'Task ',me,
     *           ' Error in doloop: stack overflow while ',
     *           'pushing wheres'
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Set the c_oploop field in the optable to indicate initialization has
c   been performed.
c---------------------------------------------------------------------------

      optable(c_oploop, iop) = 1

c---------------------------------------------------------------------------
c   Set the new values of start_op, end_op
c---------------------------------------------------------------------------

      start_op = iop
      do j = iop+1,noptable
         if (optable(c_opcode,j) .eq. enddo_op) then
               
c---------------------------------------------------------------------------
c   Does the loop index match the loop index of the DO?
c---------------------------------------------------------------------------

            call get_loop_index(optable(c_ind1,j), instr_index, dummy)
            if (instr_index .eq. index) then
               end_op = j
               if (debug) 
     *            print *,'   New start_op, end_op values = ',
     *              start_op, end_op

c----------------------------------------------------------------------------
c   Get the WHERE conditionals for the current loop.
c----------------------------------------------------------------------------

               call get_where_conditionals(optable, start_op, end_op)
               return
            endif
         endif
      enddo

      return
      end
