
c This routine adds a vector DADD to the first column of an MOIO array in
c storage. DADD must have the same number of elements as one column of the
c MOIO array.

c INPUT
c double DADD(*) : the vector of addends
c double DSCR(*) : scratch array to hold at least one column of the MOIO array
c int IDIMSCR : size of the scratch array
c int IJUNK   : (obsolete integer)
c int ILEFT   : the left  MOIO index
c int IRIGHT  : the right MOIO index

      subroutine sumsym3(dAdd,dScr,iDimScr,iJunk,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iDimScr, iJunk, iLeft, iRight
      double precision dAdd(*), dScr(iDimScr)

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer nRows, nCols, i

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert iDimScr is > 0
      if (iDimScr.lt.1) then
         print *, '@SUMSYM3: Assertion failed.'
         print *, '   iDimScr = ',iDimScr
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif

c ----------------------------------------------------------------------

      nRows = aces_list_rows(iLeft,iRight)
      nCols = aces_list_cols(iLeft,iRight)
      if ((nRows.lt.1).or.(nCols.lt.1)) return

      if (iDimScr.lt.nRows) then
         print *, '@SUMSYM3: ERROR - There is not enough scratch space',
     &            ' for one column.'
         print *, '          MOIO list = ',iLeft,iRight
         print *, '          iDimScr   = ',iDimScr
         print *, '          nRows     = ',nRows
         call aces_exit(1)
      end if

      call getlst(dScr,1,1,0,iLeft,iRight)
      do i = 1, nRows
         dScr(i) = dScr(i) + dAdd(i)
      end do
      call putlst(dScr,1,1,0,iLeft,iRight)

      return
c     end subroutine sumsym3
      end

