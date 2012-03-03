
c This routine converts seconds (in doubles) to hours, minutes, and seconds.

c INPUT
c double SECS : the number of seconds to convert

c OUTPUT
c int HH, MM : the hours and minutes in SECS
c double SS  : the remaining seconds

c#define _DEBUG_DS2HMS

      subroutine ds2hms(secs,hh,mm,ss)
      implicit none

c ARGUMENTS
      double precision secs, ss
      integer hh, mm

c PARAMETERS
      double precision s2h, s2m
      parameter (s2h=(1.d0/3600.d0), s2m=(1.d0/60.d0))

c INTERNAL VARIABLES
      double precision s
      integer h, m

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

