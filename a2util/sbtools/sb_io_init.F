
c This routine prepares the SB environment to use lists (getlst and putlst).

      subroutine sb_io_init
      implicit none

#include "aces.h"
#include "machsp.com"
#include "sbcore.com"
#include "flags.com"

      integer iMem, piocache, setptr
      logical bInitCache

      call callstack_push('SB_IO_INIT')

      bInitCache = (memknown.ne.0)
      iMem       = iintfp*iflags(37)*iflags(44)
      piocache   = setptr(1,F_INTEGER,iMem)
cTODO      call aces_io_init(iCore,piocache,iMem,bInitCache)

      call callstack_pop
      return
      end

