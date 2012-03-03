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
      subroutine handle_sp_op(op, index_table, nindex_table)
c--------------------------------------------------------------------------
c   Performs calculations on the "scratchpad".
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'scratchpad.h'
      include 'trace.h'

      integer op(loptable_entry)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer i, opcode, op1, op2, result

      integer get_symbolic_constant

      opcode = op(c_opcode)
      op1    = op(c_op1_array)
      op2    = op(c_op2_array)
      result = op(c_result_array)
      
c---------------------------------------------------------------------------
c   result field must be a valid scratchpad index.
c---------------------------------------------------------------------------

      if (result .lt. 1 .or.
     *    result .gt. mx_scratchpad) then
         print *,'Error: Invalid scratchpad address in operation.'
         print *,'Operation: ',(op(i),i=1,loptable_entry)
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Arithmetic operations:
c---------------------------------------------------------------------------

      if (opcode .eq. sp_add_op) then
         scratchpad(result) = op1 + op2
         return
      else if (opcode .eq. sp_sub_op) then
         scratchpad(result) = op1 - op2
         return
      else if (opcode .eq. sp_mult_op) then
         scratchpad(result) = op1 * op2
         return
      else if (opcode .eq. sp_div_op) then
         if (op2 .eq. 0) then
            print *,'Error: Operation contains division by 0.'
            print *,'Operation: ',(op(i),i=1,loptable_entry)
            call abort_job()
         else
            scratchpad(result) = op1 / op2
            return
         endif  
      endif
      
c---------------------------------------------------------------------
c   Logical comparison operations:
c---------------------------------------------------------------------

      if (opcode .eq. sp_equal_op) then
         if (scratchpad(op1) .eq. scratchpad(op2)) then
            scratchpad(result) = 1
         else 
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. sp_nequal_op) then
         if (scratchpad(op1) .ne. scratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. sp_ge_op) then
         if (scratchpad(op1) .ge. scratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. sp_le_op) then
         if (scratchpad(op1) .le. scratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. sp_gt_op) then
         if (scratchpad(op1) .gt. scratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      else if (opcode .eq. sp_lt_op) then
         if (scratchpad(op1) .lt. scratchpad(op2)) then
            scratchpad(result) = 1
         else
            scratchpad(result) = 0
         endif
         return
      endif
 
c--------------------------------------------------------------------------
c   Load instructions:
c--------------------------------------------------------------------------

      if (opcode .eq. sp_ldi_op) then      !     Load from instruction
         scratchpad(result) = op1
         return
      else if (opcode .eq. sp_ldi_sym_op) then ! Load from sym. table
         scratchpad(result) = get_symbolic_constant(op1)
         return
      else if (opcode .eq. sp_ldindex_op) then   ! load index segment.
         if (op1 .lt. 1 .or. 
     *       op1 .gt. nindex_table) then
            print *,'Error: Operation contains invalid index.'
            print *,'Operation: ',(op(i),i=1,loptable_entry)
            print *,'Number of index table entries = ',nindex_table
            call abort_job() 
         else
            scratchpad(result) = index_table(c_current_seg,op1)
            return
         endif
      endif
      return
      end
