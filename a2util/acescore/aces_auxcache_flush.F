
c This routine flushes any memory-resident list pointed to by quikget.
c The operation is also called a "write-back". This routine DOES NOT zero
c (destroy) the contents of quikget since doing so would perpetuate sloppy
c programming by losing the memory pointers.

      subroutine aces_auxcache_flush
      implicit none

c EXTERNAL FUNCTIONS
      integer aces_list_cols

c INTERNAL VARIABLES
      integer iFam, iGrp, iNdx, nCols

c COMMON BLOCKS
#include "icore.com" /* for icore */
#include "auxcache.com" /* for quikget */

c ----------------------------------------------------------------------

c   o write out any in-core list
      do iFam = 1, _MAX_IO_FAMS
      do iGrp = 1, _MAX_IO_GRPS
         if (quikget(iGrp,iFam).ne.0) then
            iNdx  = quikget(iGrp,iFam)
            nCols = aces_list_cols(iGrp,iFam)
            quikget(iGrp,iFam) = 0
            call putlst(icore(iNdx),1,nCols,1,iGrp,iFam)
            quikget(iGrp,iFam) = iNdx
         end if
      end do
      end do

      return
c     end subroutine aces_auxcache_flush
      end

