
c This routine changes the number of doubles per column of one array (irrep)
c in an MOIO family. The array may only have one column (like T1(A,I)).
c    This is really a shameful routine. The entire list processing subsystem
c forces the caller to retrieve at least one entire column of an array.
c However, some programmers do not feel the need to work around this
c 'feature' and instead change the size/length of the array depending
c on the data that is stored there. We use a 'private' common block
c to hold the original length of the array (set in aces_list_touch).
c If a caller attempts to resize above the original length, the program
c will die. People should consider better algorithms...

c INPUT
c int ILEFT  : the left  index in MOIO
c int IRIGHT : the right index in MOIO
c int ISIZE  : the number of doubles per column

      subroutine aces_list_resize(iLeft,iRight,iSize)
      implicit none

c ARGUMENTS
      integer iLeft, iRight, iSize

c INTERNAL VARIABLES
      integer iTmp

c COMMON BLOCKS
#include "lists.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_LIST_RESIZE: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft .lt.1).or.(_MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(_MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@ACES_LIST_RESIZE: Assertion failed.'
         print *, '   iLeft  = ',iLeft
         print *, '   iRight = ',iRight
         print *, '   MOIO rows = ',_MAX_IO_GRPS
         print *, '   MOIO cols = ',_MAX_IO_FAMS
         iTmp = 1
      end if
c   o assert iSize is whole
      if (iSize.lt.0) then
         print *, '@ACES_LIST_RESIZE: Assertion failed.'
         print *, '   iSize = ',iSize
         iTmp = 1
      end if
c   o assert the new data is <= the initialized data
      if ((iSize*moiods(iLeft,iRight)).gt.
     &         moiomxsz(iLeft,iRight)     ) then
         print *, '@ACES_LIST_RESIZE: Assertion failed.'
         print *, '   nCols = ',moiods(iLeft,iRight)
         print *, '   nRows = ',moiosz(iLeft,iRight)
         print *, '   iSize = ',iSize
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o assign the new column length
      moiosz(iLeft,iRight) = iSize

cYAU - There should be an isytyp modifier here as well; however, most one-
c      dimensional arrays are tacked onto the end of some other family.
c      isytyp(1,iRight) = ???

c   o mark MOIO as modified
      bIOMod = .true.

      return
c     end subroutine aces_list_resize
      end

