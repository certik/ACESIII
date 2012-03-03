
c This function returns the amount of memory remaining.

      integer function memleft(stack,type)
      implicit none

#include "aces.h"
#include "sb_mem.com"
#include "callstack.com"

      integer stack,type

      callstack_prev = callstack_curr
      callstack_curr='MEMLEFT'

      if (type.eq.F_INTEGER) then
         memleft=i1-iptr
      else if (type.eq.F_REAL) then
         memleft=d1-dptr
      else
         write(*,*) '@MEMLEFT: received invalid data type'
         call errex
         stop
      end if

      callstack_curr = callstack_prev
      return
      end

