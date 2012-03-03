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
      subroutine handle_sss_op(op, array_table, narray_table, 
     *                         scalar_table, nscalar_table)
c-------------------------------------------------------------------------
c   Performs all calculations of the form 
c      scalar op scalar --> scalar.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer op(loptable_entry)
      integer narray_table, nscalar_table
      integer array_table(larray_table_entry,narray_table)
      double precision scalar_table(nscalar_table)

      integer op1, op2, result
      integer iop1, iop2, iresult
      integer i, opcode, ierr
      double precision val1, val2

      op1 = op(c_op1_array)
      op2 = op(c_op2_array)
      result = op(c_result_array)
      opcode = op(c_opcode)

c--------------------------------------------------------------------------
c   Make sure all operands are scalars.
c--------------------------------------------------------------------------

      ierr = 0
      if (op2 .eq. 0) then
         if (opcode .ne. assignment_op .or.
     *       array_table(c_array_type,op1) .ne. scalar_value .or.
     *       array_table(c_array_type,result) .ne. scalar_value) 
     *      ierr = 1
      else
         if (array_table(c_array_type,op1) .ne. scalar_value .or.
     *    array_table(c_array_type,op2) .ne. scalar_value .or.
     *    array_table(c_array_type,result) .ne. scalar_value) 
     *       ierr = 1
      endif

      if (ierr .eq. 1) then
         print *,'Error in handle_sss_op: Invalid operation.'
         print *,'op = ',(op(i),i=1,loptable_entry)
         print *,'op1 type = ',array_table(c_array_type,op1)
         if (op2 .ne. 0) print *,'op2 type = ',
     *             array_table(c_array_type,op2)
         print *,'result type = ',array_table(c_array_type,result)
         call abort_job()
      endif
   
c--------------------------------------------------------------------------
c   Get operand addresses.
c--------------------------------------------------------------------------

      iop1 = array_table(c_scalar_index,op1)
      if (op2 .ne. 0) iop2 = array_table(c_scalar_index,op2)
      iresult = array_table(c_scalar_index,result)

c--------------------------------------------------------------------------
c   Perform the arithmetic.
c--------------------------------------------------------------------------

      val1 = scalar_table(iop1)
      if (op2 .ne. 0) val2 = scalar_table(iop2)

      if (opcode .eq. contraction_op) then
         scalar_table(iresult) = val1 * val2
      else if (opcode .eq. sum_op) then
         scalar_table(iresult) = val1 + val2
      else if (opcode .eq. subtract_op) then
         scalar_table(iresult) = val1 - val2
      else if (opcode .eq. assignment_op) then
         scalar_table(iresult) = val1
      else if (opcode .eq. divide_op) then
         if (val2 .eq. 0.0) then
            print *,'Error: Scalar divide by zero'
            print *,'Operation: ',(op(i),i=1,loptable_entry)
            print *,'Values = ',val1, val2
            call abort_job()
         endif

         scalar_table(iresult) = val1 / val2
      else
         print *,'Invalid scalar operation: opcode = ',opcode
         call abort_job()
      endif

      return
      end
