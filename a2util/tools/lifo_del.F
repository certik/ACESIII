
c This subroutine deletes a stack. WARNING: It is possible to delete any
c stack, but the number of active stacks will be decremented iff the id
c is that of the last stack.

c OUTPUT
c int HIGHPT : the high-water mark of used indices
c              NOTE: HighPt is not the number of touched elements unless
c              INDEX1 (from LIFO_NEW) is 1.

c INPUT
c int ID : the stack id to delete

      subroutine lifo_del(id,HighPt)
      implicit none

      integer id, HighPt

#include "lifo.com"

      if (id.lt.1.or.nStacks.lt.id) then
         print *, '@LIFO_DEL: Invalid stack id ',id
         print *, '           Must be between 1 and ',nStacks
         call c_exit(1)
      end if
      if (lifo(2,id).eq.-1) then
         print *, '@LIFO_DEL: stack ',id,' is invalid'
         call c_exit(1)
      end if

      lifo(2,id) = -1
      HighPt = lifo(3,id)-1
      if (id.eq.nStacks) nStacks = nStacks-1

      return
c     end subroutine lifo_del
      end

