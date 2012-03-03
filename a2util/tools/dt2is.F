
c This routine converts a date and time (in integers) to seconds since Jan. 1
c of THAT year, not the Epoch.

c INPUT
c int yr, mo, dy
c     hh, mm, ss : the year, month, day, hour, minute, and second

c OUTPUT
c int secs : the total seconds since 1/1

c#define _DEBUG_DT2IS

      subroutine dt2is(yr,mo,dy,hh,mm,ss,secs)
      implicit none

c ARGUMENTS
      integer yr, mo, dy, hh, mm, ss, secs

c INTERNAL VARIABLES
      integer day
      integer days(12)
      data days /0,31,59,90,120,151,181,212,243,273,304,334/

c ----------------------------------------------------------------------

#ifdef _ASSERT
c   o assert date and time are all positive
      if ((yr.lt.1).or.
     &    (mo.lt.1).or.
     &    (dy.lt.1).or.
     &    (hh.lt.0).or.
     &    (mm.lt.0).or.
     &    (ss.lt.0)    ) then
         print *, '@DT2IS: Assertion failed.'
         print *, '   yr/mo/dy = ',yr,mo,dy
         print *, '   hh:mm:ss = ',hh,mm,ss
         call c_exit(1)
      end if
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      if ((mo.gt.2).and.(iand(yr,3).eq.0)) then
         day = days(mo) + dy + 1
      else
         day = days(mo) + dy
      end if

      secs = (((day*24)+hh)*60+mm)*60+ss

      return
      end

