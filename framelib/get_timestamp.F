C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      character*24 function get_timestamp()
c-----------------------------------------------------------------------
c   Returns the current time.
c-----------------------------------------------------------------------

      implicit none

      character*8 dt
      character*10 tm

#ifdef HP
      character*24 ctime
      integer system_time
      integer time
      external time
      external ctime

      system_time = time()
      get_timestamp = ctime(system_time)
#else
      call date_and_time(dt, tm)

      get_timestamp = tm(1:2)//':'//tm(3:4)//':'//tm(5:6)
#endif
      return
      end
