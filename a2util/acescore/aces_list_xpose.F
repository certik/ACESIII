
c This routine transposes the dimensions of one MOIO array and optionally
c transposes the actual data (but not yet).

c INPUT
c int ILEFT  : the left  index in MOIO
c int IRIGHT : the right index in MOIO
c logical BXPOSE : a behavior flag
c                  = T;        transpose the data (not implemented yet)
c                  = F; do not transpose the data

      subroutine aces_list_xpose(iLeft,iRight,bXPose)
      implicit none

c ARGUMENTS
      integer iLeft, iRight
      logical bXPose

c INTERNAL VARIABLES
      integer iTmp

c COMMON BLOCKS
#include "lists.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_LIST_XPOSE: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft .lt.1).or.(_MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(_MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@ACES_LIST_XPOSE: Assertion failed.'
         print *, '   iLeft  = ',iLeft
         print *, '   iRight = ',iRight
         print *, '   MOIO rows = ',_MAX_IO_GRPS
         print *, '   MOIO cols = ',_MAX_IO_FAMS
         iTmp = 1
      end if
c   o assert bXPose is false (for now)
      if (bXPose) then
         print *, '@ACES_LIST_XPOSE: Assertion failed.'
         print *, '   bXPose = ',bXPose
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o transpose the data
      if (bXPose) then
         print *, '@ACES_LIST_XPOSE: no physical transpose yet'
      end if

c   o transpose the dimensions
      iTmp                 = moiods(iLeft,iRight)
      moiods(iLeft,iRight) = moiosz(iLeft,iRight)
      moiosz(iLeft,iRight) = iTmp

cYAU - There should be an isytyp modifier here as well; however, isytyp
c      describes an entire family and not one array.
c      isytyp(1,iRight) = ???
c      isytyp(2,iRight) = ???

c   o mark MOIO as modified
      bIOMod = .true.

      return
c     end subroutine aces_list_xpose
      end

