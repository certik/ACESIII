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
      subroutine handle_goto(optable, noptable, debug,
     *                start_op, end_op, iop)
c----------------------------------------------------------------------------
c   handle_goto: Sets the instruction pointer to the next targeted instruction
c                if the current operation is an unrestricted go_to_op.
c                Also defines a new start_op, end_op range if the target
c                instruction is not within the current range.
c
c                Also handles a jz_op instruction.  The jz_op tests the 
c                scratchpad value at the address corresponding to the 
c                result field of the instruction, and sets new start_op, 
c                end_op, and iop values if this value is 0.
c
c   Arguments:
c      optable		Operation table
c      noptable		Number of instructions in optable
c      debug		Debug flag
c      start_op         Input/output argument: see below for details.
c      end_op		Input/output argument: see below for details. 
c      iop		Input/output argument: see below for details.
c 
c   On input, start_op and end_op define the current block of insteructions
c   being executed, and iop points to the current instruction.
c
c   If the instruction opcode(c_opcode, iop) is NOT a go_to_op, start_op, 
c   end_op, and iop are not changed on output.
c
c   If the instruction IS a go_to_op, start_op and end_op may be changed to 
c   define a new loop range (if the target instruction is out of the input 
c   loop range), and iop is changed to point to the target instruction.
c-----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'scratchpad.h'
      include 'mpif.h'
      include 'parallel_info.h'

      integer noptable, start_op, end_op, iop
      integer optable(loptable_entry, noptable)
  
      integer ierr, istat, target
      integer i, op1
      integer pop_do 
      logical debug 
      
      call set_prefetch_flag(.false.)

      if (optable(c_opcode, iop) .ne. go_to_op) go to 200

c---------------------------------------------------------------------------
c   We have a go_to_op.  Find the target.
c---------------------------------------------------------------------------

      target = optable(c_result_array,iop)
      if (debug) then
         print *,'handle_goto: iop, start_op, end_op = ',
     *                 iop, start_op, end_op,' target = ',target
      endif

  100 continue
      if (target .ge. start_op .and. 
     *    target .le. end_op) then

c---------------------------------------------------------------------------
c   The jump is within the current instruction range.
c---------------------------------------------------------------------------
 
         iop = target
         if (debug) print *,'Jump to iop = ',iop,' start_op,end_op = ',
     *                       start_op,end_op
      else

c---------------------------------------------------------------------------
c   Jump is out of the current range.  Pop a new start_op, end_op
c   pair off the loop stack.
c---------------------------------------------------------------------------

         istat = pop_do(start_op, end_op)
         if (istat .lt. 0) then
            print *,'Task ',me,
     *          ' Error in handle_goto: stack underflow'
            call abort_job()
         endif

         if (debug) then
            print *,'   New start_op, end_op = ',start_op, end_op
         endif

c----------------------------------------------------------------------------
c   Re-check the target to determine if it is within the new start_op, 
c   end_op range.
c----------------------------------------------------------------------------

         go to 100
      endif

  200 continue
      if (optable(c_opcode, iop) .ne. jz_op) return

c----------------------------------------------------------------------------
c   Jump if the value is zero.
c----------------------------------------------------------------------------
         
      op1    = optable(c_op1_array,iop)
      target = optable(c_result_array, iop) 
      
      if (op1 .lt. 1 .or. op1 .gt. mx_scratchpad) then
         print *,'Error: Invalid scratchpad address in jz instruction.'
         print *,'Operation: ',(optable(i,iop),i=1,loptable_entry)
         call abort_job()
      endif

      if (target .gt. noptable) then
         print *,'Error: Invalid target address in jz instruction.'
         print *,'Operation: ',(optable(i,iop),i=1,loptable_entry)
         call abort_job()
      endif

      if (scratchpad(op1) .eq. 0) go to 100   ! reset start_op, end_op, iop
      return
      end
