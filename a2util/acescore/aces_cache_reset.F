
c This routine resets the automatic cache data. It DOES NOT write-back
c any modified slots.

      subroutine aces_cache_reset
      implicit none

c INTERNAL VARIABLES
      integer i

c COMMON BLOCKS
#include "cache.com"

c ----------------------------------------------------------------------

c   o make sure the cache subsystem is up
      if (.not.bCacheUp) return

c   o zero the cached addresses
      do i = 1, cachnum
         cachrec(i) = 0
      end do

c   o zero the LUNs
      do i = 1, cachnum
         cachfil(i) = 0
      end do

c   o reset the modification flags
      do i = 1, cachnum
         cachmod(i) = 0
      end do

c   o reset the cache event counter
      cachetime = 0

c   o zero the LRU statistics
      do i = 1, cachnum
         lrustats(i) = 0
      end do

      return
c     end subroutine aces_cache_reset
      end

