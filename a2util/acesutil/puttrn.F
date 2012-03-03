
c This routine stores the transpose of a rectangular matrix to the upper-left
c rectangular submatrix of an MOIO array. For a stored matrix
c dList(1:nRows,1:nCols) and given the input dSrc(1:iCol,1:iRow), the updated
c stored data would be dList(1:iRow:nRows,1:iCol:nCols).

c INPUT
c double DSRC(ICOL,IROW) : the transposed source array
c double DSCR(*) : scratch array to hold at least one column of the MOIO array
c int IROW   : number of rows    in the submatrix (columns in the transpose)
c int ICOL   : number of columns in the submatrix (rows    in the transpose)
c int IJUNK  : (obsolete integer)
c int ILEFT  : the left  MOIO index
c int IRIGHT : the right MOIO index

      subroutine puttrn(dSrc,dScr,iRow,iCol,iJunk,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iRow, iCol, iJunk, iLeft, iRight
      double precision dSrc(iCol,iRow), dScr(*)

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer i

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert dimensions are properly bound (zero rows or cols returns)
      if ((iRow.lt.0).or.(iRow.gt.aces_list_rows(iLeft,iRight)).or.
     &    (iCol.lt.0).or.(iCol.gt.aces_list_cols(iLeft,iRight))
     &   ) then
         print *, '@PUTTRN: Assertion failed.'
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

      if ((iRow.lt.1).or.(iCol.lt.1)) return

      if (iRow.eq.aces_list_rows(iLeft,iRight)) then
         do i = 1, iCol
            call dcopy(iRow,dSrc(i,1),iCol,dScr,1)
            call putlst(dScr,i,1,0,iLeft,iRight)
         end do
      else
         do i = 1, iCol
            call getlst(dScr,i,1,0,iLeft,iRight)
            call dcopy(iRow,dSrc(i,1),iCol,dScr,1)
            call putlst(dScr,i,1,0,iLeft,iRight)
         end do
      end if

      return
c     end subroutine puttrn
      end


