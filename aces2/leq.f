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

