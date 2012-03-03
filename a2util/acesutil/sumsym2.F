
c This routine adds an array DADD to an MOIO array in storage. DADD and the
c MOIO array must have the same number of elements.

c INPUT
c double DADD(*) : the vector of addends
c double DSCR(*) : scratch array to hold at least one column of the MOIO array
c int IDIMSCR : size of the scratch array
c int IJUNK   : (obsolete integer)
c int ILEFT   : the left  MOIO index
c int IRIGHT  : the right MOIO index

      subroutine sumsym2(dAdd,dScr,iDimScr,iJunk,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iDimScr, iJunk, iLeft, iRight
      double precision dAdd(*), dScr(iDimScr)

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer nRows, nCols, i
      integer iNdx, iTmp, iStart, nBatch, nRemain, nMax

c ----------------------------------------------------------------------

      nRows = aces_list_rows(iLeft,iRight)
      nCols = aces_list_cols(iLeft,iRight)
      if ((nRows.lt.1).or.(nCols.lt.1)) return

      if (iDimScr.lt.nRows) then
         print *, '@SUMSYM2: ERROR - There is not enough scratch space',
     &            ' for one column.'
         print *, '          MOIO list = ',iLeft,iRight
         print *, '          iDimScr   = ',iDimScr
         print *, '          nRows     = ',nRows
         call aces_exit(1)
      end if

      nRemain = nCols
      iStart  = 1
      iNdx = 1
      nMax = iDimScr/nRows
      do while (nRemain.gt.0)
         nBatch = min(nRemain,nMax)

         call getlst(dScr,iStart,nBatch,0,iLeft,iRight)
         iTmp = nRows*nBatch
         do i = 0, iTmp-1
            dScr(1+i) = dScr(1+i) + dAdd(iNdx+i)
         end do
         iNdx = iNdx + iTmp
         call putlst(dScr,iStart,nBatch,0,iLeft,iRight)

         nRemain = nRemain - nBatch
         iStart  = iStart  + nBatch
      end do

      return
c     end subroutine sumsym2
      end

