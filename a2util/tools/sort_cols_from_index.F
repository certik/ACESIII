
c This routine sorts columns according to a source vector:
c    dA(1:nRows,i) = dA(1:nRows,iFromNdx(i))
c
c On output, iFromNdx(i) will contain the same values but negated.
c Observe the iFromNdx values:
c      INPUT       OUTPUT
c    (1,2,3,4)  (-1,-2,-3,-4) : identity operation
c    (1,3,2,4)  (-1,-3,-2,-4) : 2x2 ket transpose
c    (2,3,4,1)  (-2,-3,-4,-1) : left  circular shift by one column
c    (4,1,2,3)  (-4,-1,-2,-3) : right circular shift by one column
c
c If iFromNdx(i) is zero, then the routine will not attempt to replace
c column i. There are ways to hack the inner workings using iFromNdx,
c but the safest type of sort is a 1-to-1 and onto mapping. These
c should always succede.
c
c The purpose of returning iFromNdx as a negative iToNdx vector is that
c the sort can be undone with sort_cols_to_index after negating
c the values.

c INPUT:
c int nRows = number of rows to move per column
c int nCols = number of columns/entries in dA/iFromNdx
c int lda   = number of rows in dA
c
c INPUT/OUTPUT:
c double dA(lda,nCols) = data to sort
c int iFromNdx(nCols)  = vector of source indices
c
c OUTPUT:
c double dBuf(nRows) = scratch buffer
c int iErr = return code
c            <0 -> problem with input arguments
c            =0 -> sort successful
c            >0 -> base column when an error is encountered

      subroutine sort_cols_from_index(nRows,nCols,dA,lda,dBuf,
     &                                iFromNdx,iErr)
      implicit none

c   o ARGUMENTS
      integer nRows, nCols, lda, iErr
      double precision dA(lda,nCols), dBuf(nRows)
      integer iFromNdx(nCols)

c   o VARIABLES
      integer i, pStart, iOpen
      logical bNotDone, bTrace, bVerbose
      data bTrace/.false./, bVerbose/.false./

c ----------------------------------------------------------------------

#ifdef _ASSERT
      if (lda.lt.nRows) then
         print *, '@SORT_COLS_FROM_INDEX: Assertion failed.'
         print *, '                       nRows = ',nRows
         print *, '                       lda   = ',lda
         iErr = -1
         return
      end if
#endif

      iErr = 0
      if (nRows.lt.1.or.nCols.lt.1) return

      if (bVerbose) print *, 'iFromNdx=',(iFromNdx(i),i=1,nCols)
      pStart = 1
      do while (iFromNdx(pStart).le.0.and.pStart.le.nCols)
         pStart = pStart + 1
      end do
      do while (pStart.le.nCols)
         if (bVerbose) print *, 'DEBUG: new pStart = ',pStart
         if (bVerbose) print *, 'DEBUG: ',pStart,' <- ',iFromNdx(pStart)

c      o this column must be replaced
         if (pStart.ne.iFromNdx(pStart)) then

c         o copy column pStart into buffer
            if (bTrace) print *, 'DEBUG: cp ',pStart,' buffer'
            call xcopy(nRows,dA(1,pStart),1,dBuf,1)
            iOpen = pStart
            bNotDone = .true.

            do while (bNotDone)
            if (iFromNdx(iOpen).ne.pStart) then

               if (iFromNdx(iOpen).lt.1.or.
     &             iFromNdx(iOpen).gt.nCols) then
c               o read out of bounds
                  print *, '@SORT_COLS_FROM_INDEX: read out of bounds'
                  print *, ' Column ',iOpen,' wants ',iFromNdx(iOpen),
     &                     ' out of ',nCols,' columns'
                  iErr = pStart
                  return
               end if

c            o replace the column
               if (bTrace) print *, 'DEBUG: cp ',iFromNdx(iOpen),iOpen
               call xcopy(nRows,dA(1,iFromNdx(iOpen)),1,dA(1,iOpen),1)
               iFromNdx(iOpen) = -iFromNdx(iOpen)
               iOpen = -iFromNdx(iOpen)

c           else if (iFromNdx(iOpen).eq.pStart) then
            else

               if (bTrace) print *, 'DEBUG: cp buffer',iOpen
               call xcopy(nRows,dBuf,1,dA(1,iOpen),1)
               iFromNdx(iOpen) = -iFromNdx(iOpen)
               bNotDone = .false.

c           end if (iFromNdx(iOpen).ne.pStart)
            end if
               if (bVerbose) print *,'iFromNdx=',(iFromNdx(i),i=1,nCols)
c           end do while (bNotDone)
            end do

c        end if (pStart.ne.iFromNdx(pStart))
         end if

c      o mark the replacement
         iFromNdx(pStart) = -abs(iFromNdx(pStart))
         if (bVerbose) print *, 'iFromNdx=',(iFromNdx(i),i=1,nCols)

c      o find the next column to replace
         do while (iFromNdx(pStart).le.0.and.pStart.le.nCols)
            pStart = pStart + 1
         end do

c     end do while (pStart.le.nCols)
      end do

      return
c     end subroutine sort_cols_from_index
      end

