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
      integer function str_trimlen(str)
c--------------------------------------------------------------------------
c 
c   str_trimlen: Determines the length of a FORTRAN string from the first
c                character to the last non-blank character.
c
c--------------------------------------------------------------------------
      implicit none
      character*(*) str
      integer i

      do i = len(str), 1, -1
         if (str(i:i) .ne. ' ' .and.
     *       str(i:i) .ne. char(0)) then
            str_trimlen = i
            return
         endif
      enddo

      str_trimlen = 0
      return
      end

