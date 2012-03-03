
c This subroutine creates a new LIFO indexing structure (i.e. stack).

c OUTPUT
c int ID : the new stack handle

c INPUT
c int INDEX1 : the first index to count from (cannot be -1)
c int LENGTH : the maximum number of elements in the stack
c              <  0; infinite number of elements (cf. LIFO_DEL)
c              >= 0; finite number of elements

      subroutine lifo_new(id,index1,length)
      implicit none

      integer id, index1, length

#include "lifo.com"
      external lifo_bd

      if (nStacks.eq._MAX_STACKS) then
         print *, '@LIFO_NEW: No more stacks allowed.'
         call c_exit(1)
      end if
      if (index1.eq.-1) then
         print *, '@LIFO_NEW: Invalid initial index.'
         call c_exit(1)
      end if
      if (index1.lt.-1.and.-1.le.index1+length) then
         print *, '@LIFO_NEW: Invalid index range ',index1,index1+length
         call c_exit(1)
      end if

      nStacks = nStacks+1
      id = nStacks

      lifo(1,nStacks) = index1
      lifo(2,nStacks) = index1
      lifo(3,nStacks) = index1
      lifo(4,nStacks) = index1+length-1

      return
c     end subroutine lifo_new
      end

c ----------------------------------------------------------------------

      blockdata lifo_bd
      implicit none
#undef _LIFO_COM_
#include "lifo.com"
      integer l
      parameter(l=4*_MAX_STACKS)
      data lifo /l*-1/
      data nStacks /0/
      end

