
c This routine returns the index of the first element in an integer array
c that matches the target integer.

c INPUT
c int N     : the total number of integers to test in IA
c int IA(*) : an array of integers
c int INC   : the increment (leading dimension) of IA
c int IT    : the target integer to match

      integer function iieq(n,ia,inc,it)
      implicit none

c ARGUMENTS
      integer n, ia(*), inc, it

c INTERNAL VARIABLES
      integer i, ndx

c ----------------------------------------------------------------------

      iieq = 0
      if ((n.lt.1).or.(inc.lt.1)) return

      if (inc.eq.1) then
         do i = 1, n
            if (ia(i).eq.it) then
               iieq = i
               return
            end if
         end do
      else
         ndx = 1
         do i = 1, n
            if (ia(ndx).eq.it) then
               iieq = i
               return
            end if
            ndx = ndx + inc
         end do
c     end if (inc.eq.1)
      end if

      return
      end

