
c This function releases memory that was reserved by the setptr function.
c All memory that is released is flushed with a pattern (either 0x00 or 0xFF).

      subroutine relptr(stack,type,ptr)
      implicit none

#include "aces.h"
#include "kscore.com"
#include "ks_mem.com"
#include "callstack.com"

      integer stack,type,ptr

      callstack_prev = callstack_curr
      callstack_curr = 'RELPTR'

      if (type.eq.F_INTEGER) then
         if (ptr.lt.i0 .or. ptr.gt.iptr) then
            write(*,*) '@RELPTR: invalid pointer'
            call errex
            stop
         end if
         if (memknown.ne.0) call izero(kscore(ptr),iptr-ptr)
         iptr=ptr
      else if (type.eq.F_REAL) then
         if (ptr.lt.d0 .or. ptr.gt.dptr) then
            write(*,*) '@RELPTR: invalid pointer'
            call errex
            stop
         end if
         if (memknown.ne.0) call dzero(dcore(ptr),dptr-ptr)
         dptr=ptr
      else
         write(*,*) '@RELPTR: received invalid data type'
         call errex
         stop
      end if

      callstack_curr = callstack_prev
      return
      end

