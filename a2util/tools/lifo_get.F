
c This subroutine reserves NUM elements on the stack and returns the index
c where they may be stored.

c OUTPUT
c int INDX : the stack index
c             = -1; the stack is exhausted
c            != -1; upward and onward

c INPUT
c int ID  : the stack handle (from LIFO_NEW)
c int NUM : the number of elements to reserve (0 is a valid value)

      subroutine lifo_get(id,indx,num)
      implicit none

      integer id, indx, num

#include "lifo.com"

#ifdef _ASSERT
      if (num.lt.0) then
         print *, '@LIFO_GET: Assertion failed.'
         print *, '           num = ',num
         call c_exit(1)
      end if
#endif
      if (id.lt.1.or.nStacks.lt.id) then
         print *, '@LIFO_GET: Invalid stack id ',id
         print *, '           Must be between 1 and ',nStacks
         call c_exit(1)
      end if
      if (lifo(2,id).eq.-1) then
         print *, '@LIFO_GET: stack ',id,' is invalid'
         call c_exit(1)
      end if

      indx = lifo(2,id)
      if (num.lt.1) return

      lifo(2,id) = lifo(2,id)+num
      if (lifo(1,id).le.lifo(4,id).and.lifo(4,id)+1.lt.lifo(2,id)) then
         lifo(2,id) = lifo(2,id)-num
         indx = -1
      end if
      lifo(3,id) = max(lifo(2,id),lifo(3,id))

      return
c     end subroutine lifo_get
      end

