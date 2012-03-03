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
      subroutine handle_fsp_op(op, array_table, narray_table,
     *                         scalar_table, nscalar_table_entries)
c--------------------------------------------------------------------------
c   Performs calculations on the "fscratchpad".
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'scratchpad.h'

      integer op(loptable_entry)
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer nscalar_table_entries
      double precision scalar_table(nscalar_table_entries)
      integer i, opcode, op1, op2, result, index
      double precision diff, eps

      opcode = op(c_opcode)
      op1    = op(c_op1_array)
      op2    = op(c_op2_array)
      result = op(c_result_array)
      
c---------------------------------------------------------------------------
c   result field must be a valid fscratchpad index.
c---------------------------------------------------------------------------

      if (result .lt. 1 .or.
     *    result .gt. mx_fscratchpad) then
         print *,'Error: Invalid fscratchpad address in operation.'
         print *,'Operation: ',(op(i),i=1,loptable_entry)
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Arithmetic operations:
c---------------------------------------------------------------------------

      if (opcode .eq. fl_add_op) then
         fscratchpad(result) = op1 + op2
         return
      else if (opcode .eq. fl_sub_op) then
         fscratchpad(result) = op1 - op2
         return
      else if (opcode .eq. fl_mult_op) then
         fscratchpad(result) = op1 * op2
         return
      else if (opcode .eq. fl_div_op) then
         if (op2 .eq. 0) then
            print *,'Error: Operation contains division by 0.'
            print *,'Operation: ',(op(i),i=1,loptable_entry)
            call abort_job()
         else
            fscratchpad(result) = op1 / op2
            return
         endif  
      endif
      
c---------------------------------------------------------------------
c   Logical comparison operations:
c---------------------------------------------------------------------

      eps = 1.d-10

      if (opcode .eq. fl_eq_op) then
         diff = dabs(fscratchpad(op1)-fscratchpad(op2))
         if (diff .lt. eps) then
            scratchpad(result) = 1
         else 
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. fl_ne_op) then
         diff = dabs(fscratchpad(op1)-fscratchpad(op2))
         if (diff .gt. eps) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. fl_ge_op) then
         if (fscratchpad(op1) .ge. fscratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. fl_le_op) then
         if (fscratchpad(op1) .le. fscratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. fl_gt_op) then
         if (fscratchpad(op1) .gt. fscratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. fl_lt_op) then
         if (fscratchpad(op1) .lt. fscratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      endif
 
c--------------------------------------------------------------------------
c   Load instructions:
c--------------------------------------------------------------------------

      if (opcode .eq. fl_load_value_op) then   ! load scalar value
         if (op1 .lt. 1 .or. 
     *       op1 .gt. narray_table) then
            print *,'Error: Operation contains invalid index.'
            print *,'Operation: ',(op(i),i=1,loptable_entry)
            print *,'Number of array table entries = ',narray_table
            call abort_job() 
         else
            if (array_table(c_array_type,op1) .ne. scalar_value) then
               print *,'Error: Scalar logic instruction uses a',
     *             ' non-scalar array type: array, type = ',
     *             op1, array_table(c_array_type,op1)
               call abort_job()
            endif

            index = array_table(c_scalar_index, op1)
            fscratchpad(result) = scalar_table(index)
            return
         endif
      endif
      return
      end
