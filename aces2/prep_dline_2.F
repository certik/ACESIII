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

c A debugging routine that turns dline(izl(1,i):izl(2,i)) to carets (^)
c repeating for j while making the rest spaces.

#include "iachar.h"

      subroutine prep_dline_2(dline,izl,i,j)
      implicit none

#include "linelen.par"

      character*(linelen) dline
      integer izl(2,7), i, j

      integer ndx, k

      character*1 achar, czTmp, czSpace, czCaret

c ----------------------------------------------------------------------

      czCaret = achar(_IACHAR_CARET)
      czSpace = achar(_IACHAR_SPACE)

      do ndx = 1, 7
         if (izl(1,ndx).ne.0) then
            if ((ndx.eq.i).or.(ndx.eq.j)) then
               czTmp = czCaret
            else
               czTmp = czSpace
            end if
            do k = izl(1,ndx), izl(2,ndx)
               dline(k:k) = czTmp
            end do
         end if
      end do

      return
      end

