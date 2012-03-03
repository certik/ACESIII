
c DO NOT USE THIS ROUTINE

c#define _VERBOSE_MEMSET
#define _ISCR_LEN 262144

      subroutine aces_list_memset(iLeft,iRight,iConst)
      implicit none

c ARGUMENTS
      integer iLeft, iRight, iConst

c INTERNAL VARIABLES
      integer iScr(_ISCR_LEN)
      integer iFileNum, iStat
      integer nRows, nCols
      integer iRec, iRecNdx, iTmp
      integer nLeft, nPut
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

c   o look up the dimensions of the list
      nRows = moiosz(iLeft,iRight)
      nCols = moiods(iLeft,iRight)

      if ((nRows.lt.1).or.(nCols.lt.1)) return

#ifdef _VERBOSE_MEMSET
      print *, '@ACES_LIST_MEMSET: flushing list ',iLeft,iRight,' with',
     &         iConst
#endif

c   o do a fast in-core memset to the auxiliary cache
      if (quikget(iLeft,iRight).ne.0) then
         iTmp = ifltln*nRows*nCols
         call c_memset(icore(quikget(iLeft,iRight)),iConst,iTmp)
         return
      end if

c   o memset the buffer array
      if (iprcwd.le._ISCR_LEN) then
         call c_memset(iScr,iConst,iprcln)
      else
         print *, '@ACES_LIST_MEMSET: buffer is not large enough'
         print *, '   iprcwd = ',iprcwd
         print *, '   buflen = ',_ISCR_LEN
         call aces_exit(1)
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
         print *, '@ACES_LIST_MEMSET: "',szFile(1:iLength),
     &                                   '" - file not found'
         call aces_exit(1)
      end if
#endif

c   o find the first record and integer index that point to the first element
      iRecNdx = moiowd(iLeft,iRight)
      iTmp    = (iRecNdx-1)/iprcwd
      iRec    = moio(iLeft,iRight) + iTmp
      iRecNdx = iRecNdx            - iTmp*iprcwd

c   o do the first partial record and then loop over whole records
c     which contain the remaining data
      nLeft = nRows*nCols*iintfp
      nPut  = min(nLeft,iprcwd+1-iRecNdx)
      call putlst_io(iScr,iFileNum,iRec,iRecNdx,nPut)
      nLeft = nLeft - nPut
      do while (nLeft.ne.0)
         nPut = min(nLeft,iprcwd)
         iRec = iRec + 1
         call putlst_io(iScr,iFileNum,iRec,1,nPut)
         nLeft = nLeft - nPut
      end do

#ifdef OPEN_LISTS_AS_NEEDED
c   o close the file
      close(unit=iFileNum,status='KEEP',err=666,iostat=iStat)
#endif

      return

c   o I/O error
 666  print *, '@ACES_LIST_MEMSET: I/O error'
      print *, '                   list = [',iLeft,',',iRight,']'
      print '(/)'
      call aces_io_error('ACES_LIST_MEMSET',iFileNum,iStat)

c     end subroutine aces_list_memset
      end

