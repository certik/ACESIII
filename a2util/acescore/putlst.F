
c This routine stores NCOLS columns starting at column IFIRST
c of an in-core double precision 2-dimensional array.

c INPUT
c double ISRC : the source array
c               NOTE: Internally, this is treated as an integer array
c                     since a counter runs across both ISRC and ICORE.
c int IFIRST : the index of the first column to store
c int NCOLS  : the number of columns to store
c int XCACHE : the cache buffer for the I/O operation (unused)
c              = 1; T amplitudes
c              = 2; integrals
c int ILEFT  : the left  (row) address index of the array [MOIO(ILEFT,IRIGHT)]
c int IRIGHT : the right (col) address index of the array [MOIO(ILEFT,IRIGHT)]

c#define _DEBUG_PUTLST

      subroutine putlst(iSrc,iFirst,nCols,xCache,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iSrc(*), iFirst, nCols, xCache, iLeft, iRight

c EXTERNAL FUNCTIONS
      double precision dnrm2

c INTERNAL VARIABLES
      integer iFileNum, iStat
      integer nRows
      integer iRec, iRecNdx, iTmp
      integer nLeft, nPut, iOff
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
         print *, '@PUTLST: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert iFirst > 0 and nCols >= 0
      if ((iFirst.lt.1).or.(nCols.lt.0)) then
         print *, '@PUTLST: Assertion failed.'
         print *, '   iFirst = ',iFirst
         print *, '   nCols  = ',nCols
         iTmp = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft .lt.1).or.(_MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(_MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@PUTLST: Assertion failed.'
         print *, '   iRight = ',iRight
         print *, '   iLeft  = ',iLeft
         iTmp = 1
      end if
c   o assert the list was touched
      if (moio(iLeft,iRight).lt.1) then
         print *, '@PUTLST: Assertion failed.'
         print *, '   List (',iLeft,',',iRight,') does not exist.'
         iTmp = 1
      end if
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
         print *, '@PUTLST: Error reading list (',iLeft,',',iRight,')'
         print *, '         cols available = ',iTmp
         print *, '         cols addressed = ',nCols
         print *, '         range start    = ',iFirst
         call aces_exit(1)
      end if

c   o do a fast in-core dump to the auxiliary cache
      if (quikget(iLeft,iRight).ne.0) then
         iPos = quikget(iLeft,iRight) + nRows*(iFirst-1)*iintfp
         call xcopy(nRows*nCols,iSrc,1,icore(iPos),1)
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
         print *, '@PUTLST: "',szFile(1:iLength),'" - file not found'
         call aces_exit(1)
      end if
#endif

c   o find the first record and integer index that point to the first element
      iRecNdx = moiowd(iLeft,iRight) + nRows*(iFirst-1)*iintfp
      iTmp    = (iRecNdx-1)/iprcwd
      iRec    = moio(iLeft,iRight) + iTmp
      iRecNdx = iRecNdx            - iTmp*iprcwd

#ifdef _DEBUG_PUTLST
      print *, '@PUTLST: writing cols',iFirst,' -',iFirst+nCols-1,
     &         ' of list',iLeft,',',iRight
      print *, '         file',iFileNum,', rec',iRec,', index',iRecNdx
#endif

c   o do the first partial record and then loop over whole records
c     which contain the remaining data
      nLeft = nRows*nCols*iintfp
      nPut  = min(nLeft,iprcwd+1-iRecNdx)
      call putlst_io(iSrc,iFileNum,iRec,iRecNdx,nPut)
      iOff  = 1     + nPut
      nLeft = nLeft - nPut
      do while (nLeft.ne.0)
         nPut = min(nLeft,iprcwd)
         iRec = iRec + 1
         call putlst_io(iSrc(iOff),iFileNum,iRec,1,nPut)
         iOff  = iOff  + nPut
         nLeft = nLeft - nPut
      end do

#ifdef OPEN_LISTS_AS_NEEDED
c   o close the file
      close(unit=iFileNum,status='KEEP',err=666,iostat=iStat)
#endif

      return

c   o I/O error
 666  print *, '@PUTLST: I/O error'
      print *, '         list = [',iLeft,',',iRight,']'
      print '(/)'
      call aces_io_error('PUTLST',iFileNum,iStat)

c     end subroutine putlst
      end

