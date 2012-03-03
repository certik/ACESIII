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
      subroutine nulltrm(fstr)

c------------------------------------------------------------------------
c   Null terminate a FORTRAN string for use as a string in C language
c   routines.
c
c   Arguments:
c      fstr: The string to be null-terminated.  The string must have enough
c            room for one zero byte at the end.
c-------------------------------------------------------------------------

      character*(*) fstr

c-------------------------------------------------------------------------
c   Find last non-blank character.
c-------------------------------------------------------------------------
      last = 0
      do i = len(fstr), 1, -1
         if (fstr(i:i) .ne. ' ' .and.
     *       ichar(fstr(i:i)) .ne. 0) then
            last = i
            go to 100
         endif
      enddo

  100 continue
      if (last .ne. 0 .and. last .ne. len(fstr)) then
          fstr(last+1:last+1) = char(0)
          fstr(last+2:len(fstr)) = ' '
      endif
      
      return
      end
