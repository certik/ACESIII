
c This routine converts seconds (in integers) to hours, minutes, and seconds.

c INPUT
c int SECS : the number of seconds to convert

c OUTPUT
c int HH, MM, SS : the hours, minutes, and seconds equal to SECS

c#define _DEBUG_IS2HMS

      subroutine is2hms(secs,hh,mm,ss)
      implicit none

c ARGUMENTS
      integer secs, hh, mm, ss

c PARAMETERS
      double precision s2h, s2m
      parameter (s2h=(1.d0/3600.d0), s2m=(1.d0/60.d0))

c INTERNAL VARIABLES
      integer h, m, s

c ----------------------------------------------------------------------

      s = secs
      h = s * s2h
      s = -(3600*h) + s
      m = s * s2m
      s = -(60*m) + s

      hh = h
      mm = m
      ss = s

      return
      end

