
c This subroutine sets the free-space index to INDX.

c INPUT
c int ID   : the stack handle (from LIFO_NEW)
c int INDX : the new free-space index

      subroutine lifo_set(id,indx)
      implicit none

      integer id, indx

#include "lifo.com"

      if (id.lt.1.or.nStacks.lt.id) then
         print *, '@LIFO_SET: Invalid stack id ',id
         print *, '           Must be between 1 and ',nStacks
         call c_exit(1)
      end if
      if (lifo(2,id).eq.-1) then
         print *, '@LIFO_SET: stack ',id,' is invalid'
         call c_exit(1)
      end if
      if (indx.lt.lifo(1,id).or.lifo(2,id).lt.indx) then
         print *, '@LIFO_SET: Invalid index ',indx
         print *, '           Not in [',lifo(1,id),',',lifo(2,id),']'
         call c_exit(1)
      end if

      lifo(2,id) = indx

      return
c     end subroutine lifo_set
      end

