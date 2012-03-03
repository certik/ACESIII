
c This function returns the number of doubles per column of an MOIO array.

c INPUT
c int ILEFT  : the left  index in MOIO
c int IRIGHT : the right index in MOIO

      integer function aces_list_rows(iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iLeft, iRight

c INTERNAL VARIABLES
      integer iTmp

c COMMON BLOCKS
#include "lists.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert iLeft and iRight are properly bound
      if ((iLeft .lt.1).or.(_MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(_MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@ACES_LIST_ROWS: Assertion failed.'
         print *, '   iLeft  = ',iLeft
         print *, '   iRight = ',iRight
         print *, '   MOIO rows = ',_MAX_IO_GRPS
         print *, '   MOIO cols = ',_MAX_IO_FAMS
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      aces_list_rows = moiosz(iLeft,iRight)

      return
c     end function aces_list_rows
      end

