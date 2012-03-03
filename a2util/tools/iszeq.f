
c This routine returns the index of the first element in a string array
c that matches a target string.

c NOTE: Fortran blank-pads strings and the native string comparison
c will return true if the strings are of different sizes but their
c non-blank data is equal. A different utility 'iszaleq' will return
c complete matches only.

c INPUT
c int N                : the total number of strings to test in SZA
c character*(*) SZA(*) : an array of character strings
c int INC              : the increment (leading dimension) of SZA
c character*(*) SZ     : the target string to match

      integer function iszeq(n,sza,inc,sz)
      implicit none

c ARGUMENTS
      character*(*) sza(*), sz
      integer n,inc

c INTERNAL VARIABLES
      integer i, ndx

c ----------------------------------------------------------------------

      iszeq = 0
      if ((n.lt.1).or.(inc.lt.1)) return

      if (inc.eq.1) then
         do i = 1, n
            if (sza(i).eq.sz) then
               iszeq = i
               return
            end if
         end do
      else
         ndx = 1
         do i = 1, n
            if (sza(ndx).eq.sz) then
               iszeq = i
               return
            end if
            ndx = ndx + inc
         end do
c     end if (inc.eq.1)
      end if

      return
      end

