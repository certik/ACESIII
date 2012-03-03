
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

