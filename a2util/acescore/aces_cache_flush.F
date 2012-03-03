
c This routine flushes any modified cache slots in the automatic file cache.
c The operation is also called a "write-back". This routine DOES NOT zero
c (destroy) the contents of the cache.

      subroutine aces_cache_flush
      implicit none

c INTERNAL VARIABLES
      integer i, iFile, iRec, iNdx

c COMMON BLOCKS
#include "cache.com"
#include "filspc.com" /* for iprcwd */

c ----------------------------------------------------------------------

c   o write out changed records and reset the modification flags
      if (bCacheUp) then
      do i = 1, cachnum
         if (cachmod(i).ne.0) then
            iFile = cachfil(i)
            iRec  = cachrec(i)
            iNdx  = cachndx(i)
            call aces_io_write(iFile,iRec,icache(iNdx),iprcwd)
            cachmod(i) = 0
         end if
      end do
      end if

      return
c     end subroutine aces_cache_flush
      end

