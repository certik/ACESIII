
c This routine loads a block of lists into icore-addressable memory.

      subroutine aces_auxcache_ldblk(i0,iCrSiz,
     &                               iLeft0, iLeft1, iLeftInc,
     &                               iRight0,iRight1,iRightInc)
      implicit none

c ARGUMENTS
      integer i0, iCrSiz
      integer iLeft0, iLeft1, iLeftInc, iRight0, iRight1, iRightInc

c EXTERNAL FUNCTIONS
      integer aces_list_rows, aces_list_cols

c INTERNAL VARIABLES
      integer iFam, iGrp, nRows, nCols, nInts, i

c COMMON BLOCKS
#include "icore.com" /* for icore */
#include "machsp.com" /* for iintfp */
#include "auxcache.com" /* for quikget */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert icore info is sane
      if (i0.eq.0.or.iCrSiz.lt.0) then
         print *, '@ACES_AUXCACHE_LDBLK: Assertion failed.'
         print *, '   i0     = ',i0
         print *, '   iCrSiz = ',iCrSiz
         i = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft0 .lt.1).or.(_MAX_IO_GRPS.lt.iLeft0 ).or.
     &    (iLeft1 .lt.1).or.(_MAX_IO_GRPS.lt.iLeft1 ).or.
     &    (iRight0.lt.1).or.(_MAX_IO_FAMS.lt.iRight0).or.
     &    (iRight1.lt.1).or.(_MAX_IO_FAMS.lt.iRight1).or.
     &    (iLeftInc.eq.0).or.(iRightInc.eq.0)            ) then
         print *, '@ACES_AUXCACHE_LDBLK: Assertion failed.'
         print *, '   do iFam = ',iRight0,',',iRight1,',',iRightInc
         print *, '   do iGrp = ',iLeft0 ,',',iLeft1 ,',',iLeftInc
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o load lists
      do iFam = iRight0, iRight1, iRightInc
      do iGrp = iLeft0,  iLeft1,  iLeftInc
         if (quikget(iGrp,iFam).eq.0) then
            nCols = aces_list_cols(iGrp,iFam)
            nRows = aces_list_rows(iGrp,iFam)
            nInts = iintfp*nRows*nCols
            if (0.lt.nInts.and.nInts.lt.iCrSiz) then
c               print *, '@ACES_AUXCACHE_LDBLK: Loading ',iGrp,iFam
               call getlst(icore(i0),1,nCols,1,iGrp,iFam)
               quikget(iGrp,iFam) = i0
               i0     = i0     + nInts
               iCrSiz = iCrSiz - nInts
            end if
         end if
      end do
      end do

      return
c     end subroutine aces_auxcache_ldblk
      end

