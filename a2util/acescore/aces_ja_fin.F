
c This routine finalizes the job archive environment (via the JAINDX file)
c in order for getrec and putrec to work properly for other member executables.

#include "jobarc.h"

      subroutine aces_ja_fin
      implicit none

c INTERNAL VARIABLES
      character*80 szJAINDX
      integer       iJAINDX
      integer i, iStat
      logical bOpened

c COMMON BLOCKS
#include "jobarc.com"

c ----------------------------------------------------------------------

c   o make sure the job archive subsystem is up
      if (.not.bJAUp) return

#ifdef _TAG_JA_DIRTY
c   o mark the job archive as closed
      call putrec(1,'JOBARC','JADIRTY',1,0)
#endif

#ifndef OPEN_JOBARC_AS_NEEDED
c   o close the JOBARC file
      inquire(unit=_JA_UNIT,opened=bOpened,err=666,iostat=iStat)
      if (bOpened) then
         close(unit=_JA_UNIT,status='KEEP',err=666,iostat=iStat)
      end if
#endif

c   o only update JAINDX if a record has been added
      if (bJAMod) then

c      o update JAINDX (dump the jobarc common block)
         call gfname('JAINDX',szJAINDX,iJAINDX)
         open(unit=_JA_UNIT,file=szJAINDX(1:iJAINDX),
     &        form='UNFORMATTED',status='UNKNOWN',err=666,iostat=iStat)
         rewind(_JA_UNIT,err=666,iostat=iStat)
         write(_JA_UNIT,err=666,iostat=iStat) marker, rloc, rsize, nrecs
         close(unit=_JA_UNIT,status='KEEP',err=666,iostat=iStat)

c      o reset the JOBARC modification flag
         bJAMod = .false.

c     end if (bJAMod)
      end if

#ifdef _DEBUG
c   o destroy the jobarc contents so nothing gets "lost"
      do i = 1, _MAX_JA_RECS
         marker(i) = ' '
      end do
#endif

c   o turn off the job archive subsystem flag
      bJAUp = .false.

      return

c   o I/O error
 666  print *, '@ACES_JA_FIN: I/O error'
      print '(/)'
      call aces_io_error('ACES_JA_FIN',_JA_UNIT,iStat)

c     end subroutine aces_ja_fin
      end

