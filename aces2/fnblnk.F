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

c This routine returns the index of the first non-blank character in sz.

#include "iachar.h"

       integer function fnblnk(sz)
       implicit none

       character*(*) sz
       character*1 achar, czSpace, czTab
       integer i, length

       length = len(sz)
       if (length.ne.0) then

          czSpace = achar(_IACHAR_SPACE)
          czTab   = achar(_IACHAR_TAB)

          do i = 1, length
c          o return at the first non-blank character
             if ((sz(i:i).ne.czSpace).and.
     &           (sz(i:i).ne.czTab  )     ) then
                fnblnk = i
                return
             end if
          end do

       end if

       fnblnk = 0
       return
       end

