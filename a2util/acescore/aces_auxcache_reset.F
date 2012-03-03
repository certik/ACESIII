
c This routine destroys the icore indices in quikget so aces_fin can complete.

c ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING !
c    Given the nature of the auxiliary cache, it is the programmer's
c responsibility to make sure that memory is not lost. I.e., if the memory is
c needed after the reset, then steps should be taken to record where the icore
c range starts.

      subroutine aces_auxcache_reset
      implicit none

c COMMON BLOCKS
#include "auxcache.com" /* for quikget */

c ----------------------------------------------------------------------

c   o destroy quikget
      call izero(quikget,_MAX_IO_GRPS*_MAX_IO_FAMS)

      return
c     end subroutine aces_auxcache_reset
      end

