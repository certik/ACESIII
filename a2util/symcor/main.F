
      program main
      implicit none

      integer iuhf
#include "icore.com"
#include "istart.com"

      call aces_init(icore,i0,icrsiz,iuhf,.true.)
      call symcor(icore(i0),icrsiz)
      call aces_fin

      end

