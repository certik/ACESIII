
c This routine resets the I/O statistics.

      subroutine aces_io_reset
      implicit none

c INTERNAL VARIABLES
      integer iLenMOIO, i

c COMMON BLOCKS
#include "lists.com"
#include "sympop.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_IO_RESET: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o reset the cache
      call aces_cache_reset

c   o reset the lists common block
      iLenMOIO = _MAX_IO_GRPS * _MAX_IO_FAMS
      call izero(moio,  iLenMOIO)
      call izero(moiods,iLenMOIO)
      call izero(moiosz,iLenMOIO)
      call izero(moiowd,iLenMOIO)
      call izero(moiofl,iLenMOIO)
      call izero(moiomxsz,iLenMOIO)

cYAU - There is so much bad programming related to this array that we cannot
c      destroy it.
cc   o reset the distribution types
c      call izero(isytyp,2*_MAX_IO_FAMS)

c   o reset the free-space pointers
      do i = 1, _MAX_IO_LUNS
         pRec(i)    = 1
         iIntOff(i) = 0
      end do

c   o mark MOIO as modified
      bIOMod = .true.

      return
c     end subroutine aces_io_reset
      end

