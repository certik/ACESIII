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

c This routine returns the index of the last non-blank character in sz.

#include "iachar.h"

       integer function linblnk(sz)
       implicit none

       character*(*) sz
       character*1 achar, czSpace, czTab, czNULL
       integer i, length

       length = len(sz)
       if (length.ne.0) then

          czSpace = achar(_IACHAR_SPACE)
          czTab   = achar(_IACHAR_TAB)
          czNULL  = achar(0)

          do i = length, 1, -1
c          o return at the first non-blank character
             if ((sz(i:i).ne.czSpace).and.
     &           (sz(i:i).ne.czTab  ).and.
     &           (sz(i:i).ne.czNULL )     ) then
                linblnk = i
                return
             end if
          end do

       end if

       linblnk = 0
       return
       end

