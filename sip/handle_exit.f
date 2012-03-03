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
      subroutine handle_exit(optable, noptable, debug,
     *                       array_table, narray_table,
     *                       index_table, nindex_table, 
     *                       block_map_table,
     *                       start_op, end_op, iop)
c----------------------------------------------------------------------------
c   Runtime code for the "exit" instruction: The instruction exits the 
c   innermost pardo or do loop.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'where_table.h'

      integer noptable, start_op, end_op, iop
      integer optable(loptable_entry,noptable)
      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,*)
      logical debug

      integer i, istat, opcode
      integer ierr
      integer pop_do, get_block_request
      integer result_array
      integer blkndx, get_block_number
      integer create_flag

      opcode = optable(c_opcode, end_op)
      if (opcode .ne. enddo_op .and.
     *    opcode .ne. endpardo_op) then
         print *,'Error: An exit instruction must be executed from',
     *           ' within a loop.'
         print *,'Line number is ',current_line
         print *,'start_op, end_op = ',start_op,end_op
         print *,'opcode = ',optable(c_opcode, end_op)
         call abort_job()
      endif
      
      call unset_prefetch_context()

c-------------------------------------------------------------------------
c   Reset any temp blocks created in this loop.
c-------------------------------------------------------------------------

      do i = start_op, end_op
         if (optable(c_opblock,i) .ne. 0) then

c-------------------------------------------------------------------------
c   Only clear the block_computed_flag if the array is managed by
c   blkmgr.
c-------------------------------------------------------------------------

            result_array = optable(c_result_array,i)
            if (array_table(c_array_type, result_array) .ne.
     *          static_array) then
               call block_end_of_loop(result_array,
     *               optable(c_opblock,i), optable(c_opblkndx,i),
     *               array_table, narray_table, index_table,
     *               nindex_table, block_map_table)
            endif
            optable(c_opblock,i) = 0
            optable(c_opblkndx,i) = 0
         endif
      enddo

c----------------------------------------------------------------------------
c   Pop a new set of start_op, end_op values from the loop stack.
c----------------------------------------------------------------------------

      optable(c_oploop,start_op) = 0   ! reset loop init flag
      iop = end_op + 1                 ! new iop is 1 beyond end of currentloop

c--------------------------------------------------------------------------
c   Pop a new set of iwhere, nwhere values off the stack, to properly
c   unwind the stack.
c--------------------------------------------------------------------------

            istat = pop_do(iwhere, nwhere)
            if (istat .lt. 0) then
               print *,'Task ',me,
     *          ' Error in handle_exit: stack underflow'
               call abort_job()
            endif

      istat = pop_do(start_op, end_op)
      if (istat .lt. 0) then
         print *,'Task ',me,
     *        ' Error in handle_exit: stack underflow'
         call abort_job()
      endif

      if (start_op .gt. noptable .or.
     *       start_op .lt. 0 .or.
     *       end_op .gt. noptable .or.
     *       end_op .lt. 0) then
         print *,'Task ',me,' handle_exit: start_op, end_op ',
     *       start_op, end_op,' noptable ',noptable
         print *,'New iop = ',iop
         call abort_job()
      endif

      return
      end
