
c This routine returns .TRUE. if the two strings are NOT lexically equal.
c WARNING: Since Fortran pads the end of all strings with spaces,
c be sure to pass only the relevant substring.

      logical function lne(s1,s2)
      implicit none
      character*(*) s1, s2
      logical leq
      lne = .not.leq(s1,s2)
      return
      end

