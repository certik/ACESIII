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
      subroutine handle_call(optable, noptable, proctab, debug,
     *                start_op, end_op, iop)
c----------------------------------------------------------------------------
c   handle_call: Saves the current start_op, end_op and iop on the stack,
c                then determines a new start_op, end_op corresponding to 
c                target proc.  Finally, the new value of iop is set to 
c                the first instruction of the target proc.
c
c   Arguments:
c      optable		Operation table
c      noptable		Number of instructions in optable
c      debug		Debug flag
c      start_op         Input/output argument
c      end_op		Input/output argument
c      iop		Input/output argument
c 
c-----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer noptable, start_op, end_op, iop
      integer optable(loptable_entry, noptable)
      integer proctab(2,*)
      integer timer
  
      integer i, ierr, istat, target
      integer push_do 
      logical debug 
      
      if (optable(c_opcode, iop) .ne. call_op) return

c----------------------------------------------------------------------------
c   Start the procedure timer.
c----------------------------------------------------------------------------

      timer = optable(c_op1_array,iop)
      call timer_start(timer)

c---------------------------------------------------------------------------
c   We have a call_op.  Find the target.
c---------------------------------------------------------------------------

      target = optable(c_result_array,iop)
      if (trace .and. and(tracelevel, proc_trace) .ne. 0) then
         print *,'Task ',me,' Calling procedure at line ',
     *       optable(c_lineno,target),' from line ',
     *       optable(c_lineno,iop) 
      endif

c----------------------------------------------------------------------------
c   Save return values on stack.
c----------------------------------------------------------------------------

      istat = push_do(start_op, end_op)
      if (istat .lt. 0) then
         print *,'Task ',me,
     *       ' Error in handle_call: stack overflow'
         call abort_job()
      endif

      istat = push_do(iop+1, timer)
      if (istat .lt. 0) then
         print *,'Task ',me,
     *       ' Error in handle_call: stack overflow'
         call abort_job()
      endif

c----------------------------------------------------------------------------
c   Determine the new program range.
c----------------------------------------------------------------------------

      start_op = target
      end_op   = proctab(2,1)
      iop      = start_op

      return
      end
