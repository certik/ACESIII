
c This routine returns .TRUE. if the two strings are lexically equal.
c WARNING: Since Fortran pads the end of all strings with spaces,
c be sure to pass only the relevant substring.

      logical function leq(s1,s2)
      implicit none

      character*(*) s1, s2

      integer length, c_strncasecmp

      length = len(s1)
      if (length.ne.len(s2)) then
         leq = .false.
      else
         leq = (c_strncasecmp(s1,s2,length).eq.0)
      end if

      return
      end

