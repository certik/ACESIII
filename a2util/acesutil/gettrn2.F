
c This routine retrieves the transpose of a lower-left rectangular submatrix
c of an MOIO array. For a stored matrix dList(1:nRows,1:nCols) and given
c boundaries iRow and iCol, the output would be dDest(1:iCol,iRow:nRows).

c INPUT
c double DSCR(*) : scratch array to hold at least one column of the MOIO array
c int IROW   : first row    in the submatrix (column in the transpose)
c int ICOL   : last  column in the submatrix (row    in the transpose)
c int IJUNK  : (obsolete integer)
c int ILEFT  : the left  MOIO index
c int IRIGHT : the right MOIO index

c OUTPUT
c double DDEST(*) : the destination array for the transpose

      subroutine gettrn2(dDest,dScr,iRow,iCol,iJunk,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iRow, iCol, iJunk, iLeft, iRight
      double precision dDest(iCol,*), dScr(*)

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer nRows, i, iTmp

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert dimensions are properly bound (zero cols returns)
      if ((iRow.lt.1).or.(iRow.gt.aces_list_rows(iLeft,iRight)).or.
     &    (iCol.lt.0).or.(iCol.gt.aces_list_cols(iLeft,iRight))
     &   ) then
         print *, '@GETTRN2: Assertion failed.'
         print *, '   MOIO list = ',iLeft,iRight
         print *, '   MOIO rows = ',aces_list_rows(iLeft,iRight)
         print *, '   MOIO cols = ',aces_list_cols(iLeft,iRight)
         print *, '   iRow = ',iRow
         print *, '   iCol = ',iCol
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif

c ----------------------------------------------------------------------

      nRows = aces_list_rows(iLeft,iRight)
      if ((iRow.gt.nRows).or.(iCol.lt.1)) return
      if (iRow.lt.1) then
         print *, '@GETTRN2: ERROR - the valid row domain is [ 1,',
     &                                                         nRows,']'
         print *, '   MOIO list = ',iLeft,iRight
         print *, '   MOIO rows = ',aces_list_rows(iLeft,iRight)
         print *, '   MOIO cols = ',aces_list_cols(iLeft,iRight)
         print *, '   iRow = ',iRow
         print *, '   iCol = ',iCol
         call aces_exit(1)
      end if

      iTmp = nRows+1-iRow
      do i = 1, iCol
         call getlst(dScr,i,1,0,iLeft,iRight)
         call dcopy(iTmp,dScr(iRow),1,dDest(i,1),iCol)
      end do

      return
c     end subroutine gettrn
      end

