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

cSG/JDW 6/6/95
c Change of name for compatibility with newer IBM operating systems (XLF 3.2).

      character*(*) function trmblk(sz)
      implicit none

      character*(*) sz
      integer i, fnblnk, linblnk, iStart, iEnd

      trmblk=' '

      iStart = fnblnk(sz)
      if (iStart.ne.0) then
         iEnd = min(linblnk(sz),iStart-1+len(trmblk))
         trmblk = sz(iStart:iEnd)
      end if

      return
      end

