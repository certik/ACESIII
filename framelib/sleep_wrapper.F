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
      subroutine sleep_wrapper(nsecs)
c---------------------------------------------------------------------------
c   Wrapper routine for machine-dependent sleep routines.
c---------------------------------------------------------------------------
      implicit none
      integer nsecs
      integer*4 i4sleep,i4junk,i4err
#ifdef __crayx1
      i4sleep = nsecs
      call pxfsleep(i4sleep,i4junk,i4err)
#else
      call sleep(nsecs)
#endif
      return
      end
