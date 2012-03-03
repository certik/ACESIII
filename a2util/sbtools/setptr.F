
c This function returns a pointer to some location in an array (either
c kscore or dcore).  It also keeps track of memory usage and memory
c requirements.
c
c The first argument (stack) is for future use.  Multiple stacks will
c be implemented for multithread processing.

      integer function setptr(stack,type,num)
      implicit none

#include "aces.h"
#include "machsp.com"
#include "ks_mem.com"
#include "kscore.com"
#include "callstack.com"

      integer stack,type,num

      callstack_prev = callstack_curr
      callstack_curr = 'SETPTR'

      if (type.eq.F_INTEGER) then
         setptr=iptr
         iptr=iptr+num+iand(num,iintfp-1)
      else if (type.eq.F_REAL) then
         setptr=dptr
         dptr=dptr+num
      else
         write(*,*) '@SETPTR: received invalid data type'
         call errex
         stop
      end if

      if (memknown.eq.0) then
         ineeded=max(ineeded,iptr-1)
         dneeded=max(dneeded,dptr-1)
      else
         if (iptr.gt.i1 .or. dptr.gt.d1) then
            write(*,*) '@SETPTR: ran out of memory'
            call errex
            stop
         end if
      end if

      callstack_curr = callstack_prev
      return
      end

