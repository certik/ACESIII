
c This routine retrieves NCOLS columns starting at column IFIRST
c of a stored double precision 2-dimensional array.

c OUTPUT
c double IDEST : the destination array
c                NOTE: Internally, this is treated as an integer array
c                      since a counter runs across both IDEST and ICORE.

c INPUT
c int IFIRST : the index of the first column to retrieve
c int NCOLS  : the number of columns to retrieve
c int XCACHE : (OBSOLETE) the cache buffer for the I/O operation
c              = 1; T amplitudes
c              = 2; integrals
c int ILEFT  : the left  (row) address index of the array [MOIO(ILEFT,IRIGHT)]
c int IRIGHT : the right (col) address index of the array [MOIO(ILEFT,IRIGHT)]

      subroutine getlst(iDest,iFirst,nCols,xCache,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iDest(*), iFirst, nCols, xCache, iLeft, iRight

c EXTERNAL FUNCTIONS
      double precision dnrm2, dTmp
      integer idamax

c INTERNAL VARIABLES
      integer iFileNum, iStat
      integer nRows
      integer iRec, iRecNdx, iTmp
      integer nLeft, nGet, iOff
      integer iPos
#ifdef OPEN_LISTS_AS_NEEDED
#include "szfiles.dat"
      integer ihFile
      character*80 szFile
      integer      iLength
      logical bExist
#endif

c COMMON BLOCKS
#include "icore.com" /* for icore */
#include "lists.com" /* for moio info */
#include "machsp.com" /* for iintfp */
#include "filspc.com" /* for iprcwd */
#include "auxcache.com" /* for quikget */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@GETLST: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert iFirst > 0 and nCols >= 0
      if ((iFirst.lt.1).or.(nCols.lt.0)) then
         print *, '@GETLST: Assertion failed.'
         print *, '   iFirst = ',iFirst
         print *, '   nCols  = ',nCols
         iTmp = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft .lt.1).or.(_MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(_MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@GETLST: Assertion failed.'
         print *, '   iRight = ',iRight
         print *, '   iLeft  = ',iLeft
         iTmp = 1
      end if
c   o assert the list was touched
      if (moio(iLeft,iRight).lt.1) then
         print *, '@GETLST: Assertion failed.'
         print *, '   List (',iLeft,',',iRight,') does not exist.'
         iTmp = 1
      end if
c   o temporary trap to gauge list reads
c
c      if (iRight .EQ. 17 .or. iRight .eq. 18 .or.
c     &    iRight .EQ. 19 .or. iRight .eq. 20 .or.
c     &    iRight .EQ. 21 .or. iRight .eq. 22) then
c           Print *, 'Reading Redundent lists'
c           iTmp = 1
c      endif 
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o look up the length of one column
      nRows = moiosz(iLeft,iRight)

      if ((nRows.lt.1).or.(nCols.lt.1)) return

c   o make sure the column range requested is properly bound
      iTmp = moiods(iLeft,iRight)
      if ((iFirst.lt.1).or.(iTmp.lt.iFirst).or.
     &    (iTmp.lt.(iFirst-1+nCols))           ) then
         print *, '@GETLST: Error reading list (',iLeft,',',iRight,')'
         print *, '         cols available = ',iTmp
         print *, '         cols requested = ',nCols
         print *, '         range start    = ',iFirst
         call aces_exit(1)
      end if

c   o do a fast in-core fetch from the auxiliary cache
      if (quikget(iLeft,iRight).ne.0) then
         iPos = quikget(iLeft,iRight) + nRows*(iFirst-1)*iintfp
         call xcopy(nRows*nCols,icore(iPos),1,iDest,1)
         return
      end if

c   o determine which file this array lives in
      iFileNum = moiofl(iLeft,iRight)

#ifdef OPEN_LISTS_AS_NEEDED
c   o open the file that contains this array
      ihFile = 1+iFileNum-_FIRST_IO_LUN
      call gfname(szFiles(ihFile),szFile,iLength)
      inquire(file=szFile(1:iLength),exist=bExist,err=666,iostat=iStat)
      if (bExist) then
         call aces_io_open(ihFile)
      else
         print *, '@GETLST: "',szFile(1:iLength),'" - file not found'
         call aces_exit(1)
      end if
#endif

c   o find the first record and integer index that point to the first element
      iRecNdx = moiowd(iLeft,iRight) + nRows*(iFirst-1)*iintfp
      iTmp    = (iRecNdx-1)/iprcwd
      iRec    = moio(iLeft,iRight) + iTmp
      iRecNdx = iRecNdx            - iTmp*iprcwd

c   o do the first partial record and then loop over whole records
c     which contain the remaining data
      nLeft = nRows*nCols*iintfp
      nGet  = min(nLeft,iprcwd+1-iRecNdx)
      call getlst_io(iDest,iFileNum,iRec,iRecNdx,nGet)
      iOff  = 1     + nGet
      nLeft = nLeft - nGet
      do while (nLeft.ne.0)
         nGet = min(nLeft,iprcwd)
         iRec = iRec + 1
         call getlst_io(iDest(iOff),iFileNum,iRec,1,nGet)
         iOff  = iOff  + nGet
         nLeft = nLeft - nGet
      end do

#ifdef OPEN_LISTS_AS_NEEDED
c   o close the file
      close(unit=iFileNum,status='KEEP',err=666,iostat=iStat)
#endif

#ifdef _DEBUG_ACES_IO
c    o trap for NaN
       do iTmp = 1, nRows*nCols*iintfp, iintfp
          call nan_abort(iDest(iTmp))
       end do
c       print *, '@GETLST: (DEBUG) inverse norm of list ',iLeft,iRight
c       print *, '         1/dnrm2 = ',1.d0/dnrm2(nRows*nCols,iDest,1)
c       print *, '@GETLST: (DEBUG) norm of list ',iLeft,iRight
c       print *, '         dnrm2 = ',dnrm2(nRows*nCols,iDest,1)
c       print *, '@GETLST: (DEBUG) max of list ',iLeft,iRight
c       iTmp = idamax(nRows*nCols,iDest,1)
c       dTmp = iDest(iintfp*(iTmp-1)+1)
c       call nan_abort(dTmp)
c      print *, '@GETLST: (args)',iFirst,nCols,' (',iLeft,iRight,')'
c      print *, '         (dims)',nRows,moiods(iLeft,iRight)
c      print *, '         (addr)',moio(iLeft,iRight),moiowd(iLeft,iRight)
c      call dmat_prt(iDest,nRows,nCols,nRows)
c      print *, '         (norm)',dnrm2(nRows*nCols,iDest,1)
#endif

      return

c   o I/O error
 666  print *, '@GETLST: I/O error'
      print *, '         list = [',iLeft,',',iRight,']'
      print '(/)'
      call aces_io_error('GETLST',iFileNum,iStat)

c     end subroutine getlst
      end

