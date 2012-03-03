
c This routine is a wrapper for looping over sumsym2. (I know, don't ask,
c it's a legacy thing...) It adds a group of arrays DADD to an MOIO family in
c storage. DADD and the MOIO family (up to the number of irreps) must have the
c same number of elements.

c INPUT
c double DADD(*) : the vector of addends
c double DSCR(*) : scratch array to hold at least one column of any MOIO array
c int IDIMSCR : size of the scratch array
c int IFAMILY : the right MOIO index

      subroutine sumsym(dAdd,dScr,iDimScr,iFamily)
      implicit none

c ARGUMENTS
      integer iDimScr, iFamily
      double precision dAdd(*), dScr(iDimScr)

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer iGroup, iNdx, nRows, nCols

c COMMON BLOCKS
#include "syminf.com" /* for nirrep */

c ----------------------------------------------------------------------

      iNdx = 1
      do iGroup = 1, nirrep
         call sumsym2(dAdd(iNdx),dScr,iDimScr,0,iGroup,iFamily)
         nRows = aces_list_rows(iGroup,iFamily)
         nCols = aces_list_cols(iGroup,iFamily)
         iNdx = iNdx + (nRows*nCols)
      end do

      return
c     end subroutine sumsym
      end

