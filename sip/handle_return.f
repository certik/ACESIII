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
      subroutine handle_return(optable, noptable, debug,
     *                start_op, end_op, iop)
c----------------------------------------------------------------------------
c   handle_return: Replaces start_op, end_op and iop with values pulled 
c                  from the stack.  This is equivalent to a procedure
c                  return.
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
  
      integer ierr, istat, timer
      integer pop_do 
      logical debug 
      
      if (optable(c_opcode, iop) .ne. return_op) return

      if (trace .and. and(tracelevel, proc_trace) .ne. 0) then
         print *,'Task ',me,' Return from proc at line ',
     *            optable(c_lineno,iop)
      endif

c---------------------------------------------------------------------------
c   We have a return_op.  Find the return address.
c---------------------------------------------------------------------------

      istat = pop_do(iop, timer)
      if (istat .lt. 0) then
         print *,'Task ',me,
     *      ' Error in handle_return: stack underflow'
         call abort_job()
      endif

      istat = pop_do(start_op, end_op)
      if (istat .lt. 0) then
         print *,'Task ',me,
     *      ' Error in handle_return: stack underflow'
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Stop the procedure timer.
c---------------------------------------------------------------------------

      if (timer .lt. 1) then
         print *,'Task ',me,' Timer for iop ',iop,' is ',timer
         call abort_job()
      endif  
      call update_timer(timer)

      return
      end
