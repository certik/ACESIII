
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

