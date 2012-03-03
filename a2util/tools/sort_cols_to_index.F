
c This routine sorts columns according to a destination vector:
c    dA(1:nRows,iToNdx(i)) = dA(1:nRows,i)
c
c On output, iToNdx(i) will contain the negated index of the presorted
c column that now occupies column i. Observe the iToNdx values:
c      INPUT       OUTPUT
c    (1,2,3,4)  (-1,-2,-3,-4) : identity operation
c    (1,3,2,4)  (-1,-3,-2,-4) : 2x2 ket transpose
c    (2,3,4,1)  (-4,-1,-2,-3) : right circular shift by one column
c    (4,1,2,3)  (-2,-3,-4,-1) : left  circular shift by one column
c
c If iToNdx(i) is zero, then the routine will not attempt to move column
c i and it will allow another column to overwrite it (according to
c iToNdx, of course). There are ways to hack the inner workings using
c iToNdx, but the safest type of sort is a 1-to-1 and onto mapping.
c These should always succede.
c
c The purpose of returning iToNdx as a negative iToNdx vector is that
c the sort can be undone with sort_cols_to_index after negating
c the values.

c INPUT:
c int nRows = number of rows to move per column
c int nCols = number of columns/entries in dA/iToNdx
c int lda   = number of rows in dA
c
c INPUT/OUTPUT:
c double dA(lda,nCols) = data to sort
c int iToNdx(nCols)    = vector of destination indices
c
c OUTPUT:
c double dBuf(nRows) = scratch buffer
c int iErr = return code
c            <0 -> problem with input arguments
c            =0 -> sort successful
c            >0 -> base column when an error is encountered

      subroutine sort_cols_to_index(nRows,nCols,dA,lda,dBuf,iToNdx,iErr)
      implicit none

c   o ARGUMENTS
      integer nRows, nCols, lda, iErr
      double precision dA(lda,nCols), dBuf(nRows)
      integer iToNdx(nCols)

c   o VARIABLES
      integer i, pStart, iBufTo, iOpen, iFrom
      logical bNotDone, bTrace, bVerbose
      data bTrace/.false./, bVerbose/.false./

c ----------------------------------------------------------------------

#ifdef _ASSERT
      if (lda.lt.nRows) then
         print *, '@SORT_COLS_TO_INDEX: Assertion failed.'
         print *, '                     nRows = ',nRows
         print *, '                     lda   = ',lda
         iErr = -1
         return
      end if
#endif

      iErr = 0
      if (nRows.lt.1.or.nCols.lt.1) return

      if (bVerbose) print *, 'iToNdx=',(iToNdx(i),i=1,nCols)
      pStart = 1
      do while (iToNdx(pStart).le.0.and.pStart.le.nCols)
         pStart = pStart + 1
      end do
      do while (pStart.le.nCols)
         if (bVerbose) print *, 'DEBUG: new pStart = ',pStart
         if (bVerbose) print *, 'DEBUG: ',pStart,' -> ',iToNdx(pStart)

         iBufTo = iToNdx(pStart)

c      o this column must be moved
         if (pStart.ne.iToNdx(pStart)) then

            if (iToNdx(pStart).gt.nCols) then
c            o write out of bounds
               print *, '@SORT_COLS_TO_INDEX: write out of bounds'
               print *, ' Column ',pStart,' sent to ',iToNdx(pStart),
     &                  ' out of ',nCols,' columns'
               iErr = pStart
               return
            end if

            if (iToNdx(iToNdx(pStart)).gt.0) then
c            o copy column pStart into buffer
               if (bTrace) print *, 'DEBUG: cp ',pStart,' buffer'
               call xcopy(nRows,dA(1,pStart),1,dBuf,1)
               iOpen = pStart
               iToNdx(pStart) = 0
               bNotDone = .true.
            else
               if (iToNdx(iToNdx(pStart)).eq.0) then
c               o move immediately
                  iOpen = iToNdx(pStart)
                  if (bTrace) print *, 'DEBUG: cp ',pStart,iOpen
                  call xcopy(nRows,dA(1,pStart),1,dA(1,iOpen),1)
                  iToNdx(pStart) = 0
                  bNotDone = .false.
               else
c               o pStart is replacing a column that has moved already
                  print *, '@SORT_COLS_TO_INDEX: Output dependence.'
                  print *, ' Collision at column ',iToNdx(pStart),
     &                     ' from columns ',-iToNdx(iToNdx(pStart)),
     &                     ' and ',pStart
                  iErr = pStart
                  return
               end if
            end if
            if (bVerbose) print *, 'iToNdx=',(iToNdx(i),i=1,nCols)

            do while (bNotDone)
            if (iOpen.ne.iBufTo) then

c            o find the column that belongs in the open column
               iFrom = pStart + 1
               do while (iToNdx(iFrom).ne.iOpen.and.iFrom.le.nCols)
                  iFrom = iFrom + 1
               end do
               if (iFrom.gt.nCols) then
                  print *, '@SORT_COLS_TO_INDEX: buffer occupied'
                  print *, ' No column is scheduled to replace column ',
     &                     iOpen
                  iErr = pStart
                  return
               end if

c            o move the column
               if (bTrace) print *, 'DEBUG: cp ',iFrom,iOpen
               call xcopy(nRows,dA(1,iFrom),1,dA(1,iOpen),1)
               iToNdx(iOpen) = -iFrom
               iToNdx(iFrom) = 0
               iOpen = iFrom
               if (bVerbose) print *, 'iToNdx=',(iToNdx(i),i=1,nCols)

c           else if (iOpen.eq.iBufTo) then
            else

               if (bTrace) print *, 'DEBUG: cp buffer',iOpen
               call xcopy(nRows,dBuf,1,dA(1,iOpen),1)
               bNotDone = .false.

c           end if (iOpen.ne.iBufTo)
            end if
c           end do while (bNotDone)
            end do

c        end if (pStart.ne.iToNdx(pStart))
         end if

c      o mark the move
         iToNdx(iBufTo) = -pStart
         if (bVerbose) print *, 'iToNdx=',(iToNdx(i),i=1,nCols)

c      o find the next column to move
         do while (iToNdx(pStart).le.0.and.pStart.le.nCols)
            pStart = pStart + 1
         end do

c     end do while (pStart.le.nCols)
      end do

      return
c     end subroutine sort_cols_to_index
      end

