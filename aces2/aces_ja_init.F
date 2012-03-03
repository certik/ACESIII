C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.

c This routine initializes the job archive environment (via the JAINDX file)
c in order for getrec and putrec to work properly.

#include "jobarc.h"

      subroutine aces_ja_init
      implicit none

c INTERNAL VARIABLES
      character*80 szJOBARC, szJAINDX
      integer       iJOBARC,  iJAINDX
      integer i, iStat
      logical bExist, bOpened
      integer*8 iBuf(_JA_RECWD)

c COMMON BLOCKS
#include "jobarc.com"
#include "icdacc.com" /* for idaccm */

c ----------------------------------------------------------------------

c   o die instead of silently return if the job archive subsystem is already up
      if (bJAUp) then
         print *, '@ACES_JA_INIT: The job archive subsystem is already',
     &            ' initialized.'
c         call aces_exit(1)
         call abort_job()
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
            call abort_job() 
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
c         call izero(iBuf,_JA_RECWD)
         do i = 1, _JA_RECWD
            iBuf(i) = 0
         enddo

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
      call igetrec(-1,'JOBARC','JADIRTY',1,i)
      if (i.ne.0) then
         print '(/)'
         print *, '@ACES_JA_INIT: WARNING - The job archive was not ',
     &            'finalized by the previous'
         print *, '               ACES Member Executable. Any records ',
     &            'added by that process have'
         print *, '               been lost.'
         print '(/)'
      else
         call iputrec(1,'JOBARC','JADIRTY',1,1)
      end if
#endif

      return

c   o I/O error
 666  print *, '@ACES_JA_INIT: I/O error ',iStat
      call abort_job()

c     end subroutine aces_ja_init
      end

