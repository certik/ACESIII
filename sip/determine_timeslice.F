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
      subroutine determine_timeslice(timeslice, timer_count)
c--------------------------------------------------------------------------
c   Returns system-dependent vallues for timeslice and timer_count 
c   arguments.
c--------------------------------------------------------------------------
      implicit none
      double precision timeslice
      integer timer_count
      integer itimeslice

#ifdef __crayx1
      itimeslice = 1000   ! 1000 microseconds
#else
      itimeslice = 500
#endif

#ifdef ALTIX
      timer_count = 1
#else
      timer_count = 10
#endif
      timeslice = itimeslice * 1.e-6   ! convert microseconds to seconds.
      return
      end

