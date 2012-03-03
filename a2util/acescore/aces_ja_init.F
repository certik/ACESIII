
c This routine initializes the job archive environment (via the JAINDX file)
c in order for getrec and putrec to work properly.

#include "jobarc.h"

      subroutine aces_ja_init
      implicit none

c INTERNAL VARIABLES
      character*80 szJOBARC, szJAINDX
      integer       iJOBARC,  iJAINDX
      integer iBuf(_JA_RECWD), i, iStat
      logical bExist, bOpened

c COMMON BLOCKS
#include "jobarc.com"
#include "icdacc.com" /* for idaccm */

c ----------------------------------------------------------------------

c   o die instead of silently return if the job archive subsystem is already up
      if (bJAUp) then
         print *, '@ACES_JA_INIT: The job archive subsystem is already',
     &            ' initialized.'
         call aces_exit(1)
      end if

c   o turn on the job archive subsystem flag
      bJAUp = .true.

c   o get the external JOBARC file name
      call gfname('JOBARC',szJOBARC,iJOBARC)

c   o initialize the jobarc common block
      call gfname('JAINDX',szJAINDX,iJAINDX)
      inquire(file=szJAINDX(1:iJAINDX),exist=bExist,
     &        err=666,iostat=iStat)
      if (bExist) then
c      o JOBARC had better exist
         inquire(file=szJOBARC(1:iJOBARC),exist=bExist,
     &           err=666,iostat=iStat)
         if (.not.bExist) then
            print *, '@ACES_JA_INIT: JOBARC does not exist'
            call aces_exit(1)
         end if
c      o process JAINDX
         open(unit=_JA_UNIT,file=szJAINDX(1:iJAINDX),
     &        form='UNFORMATTED',status='OLD',err=666,iostat=iStat)
         rewind(_JA_UNIT,err=666,iostat=iStat)
         read(_JA_UNIT,err=666,iostat=iStat) marker, rloc, rsize, nrecs
         close(unit=_JA_UNIT,status='KEEP',err=666,iostat=iStat)
      else
c      o reset the records
         do i = 1, _MAX_JA_RECS
            marker(i) = 'OPENSLOT'
         end do
         rloc(1) = 1
         nrecs   = 0
      end if
      irecwd = _JA_RECWD
      irecln = _JA_RECWD*idaccm

c   o condition the JOBARC file
      inquire(file=szJOBARC(1:iJOBARC),exist=bExist,opened=bOpened,
     &        err=666,iostat=iStat)
      if (bExist) then
#ifndef OPEN_JOBARC_AS_NEEDED
c      o open the JOBARC file
         if (.not.bOpened) then
            open(unit=_JA_UNIT,file=szJOBARC(1:iJOBARC),
     &           form='UNFORMATTED',access='DIRECT',
     &           status='OLD',recl=irecln,err=666,iostat=iStat)
         end if
#endif
      else
c      o write out one all-zero record (bug fix: Ajith 03/25/97)
         call izero(iBuf,_JA_RECWD)
         open(unit=_JA_UNIT,file=szJOBARC(1:iJOBARC),
     &        form='UNFORMATTED',access='DIRECT',
     &        status='NEW',recl=irecln,err=666,iostat=iStat)
         write(unit=_JA_UNIT,rec=1,err=666,iostat=iStat) iBuf
#ifdef OPEN_JOBARC_AS_NEEDED
         close(_JA_UNIT,err=666,iostat=iStat)
#endif
      end if

#ifdef _TAG_JA_DIRTY
c   o check if the job archive was properly closed and mark it as open
      call getrec(-1,'JOBARC','JADIRTY',1,i)
      if (i.ne.0) then
         print '(/)'
         print *, '@ACES_JA_INIT: WARNING - The job archive was not ',
     &            'finalized by the previous'
         print *, '               ACES Member Executable. Any records ',
     &            'added by that process have'
         print *, '               been lost.'
         print '(/)'
      else
         call putrec(1,'JOBARC','JADIRTY',1,1)
      end if
#endif

      return

c   o I/O error
 666  print *, '@ACES_JA_INIT: I/O error'
      print '(/)'
      call aces_io_error('ACES_JA_INIT',_JA_UNIT,iStat)

c     end subroutine aces_ja_init
      end

