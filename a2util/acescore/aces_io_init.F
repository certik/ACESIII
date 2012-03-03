
c This routine initializes the I/O subsystem.

c INPUT
c int ICORE : the core memory anchor
c int I0    : the index of the first usable core integer
c int IMEM  : the integer-length of the core (passed to aces_cache_init)
c logical BINITCACHE : do (not) initialize the I/O cache subsystem at iCore(i0)

      subroutine aces_io_init(iCore,i0,iMem,bInitCache)
      implicit none

c ARGUMENTS
      integer iCore(*), i0, iMem
      logical bInitCache

c EXTERNAL FUNCTIONS
      integer iieq

c INTERNAL VARIABLES
#include "szfiles.dat"
      character*80 szFile
      integer      iLength
      integer iLenMOIO, ihFile, iUnit, iPos, iTmp, iStat
      integer nprcwd, nslots
      logical bExist

c COMMON BLOCKS
#include "lists.com"
#include "sympop.com" /* for isytyp */
#include "auxcache.com"
#include "flags.com" /* for iflags(37) and iflags(44) */
#include "machsp.com" /* for ibitwd */

c ----------------------------------------------------------------------

c   o die instead of silently return if the I/O subsystem is already up
      if (bIOUp) then
         print *, '@ACES_IO_INIT: The I/O subsystem is already',
     &            ' initialized.'
         call aces_exit(1)
      end if

c   o turn on the I/O subsystem flag
      bIOUp = .true.

c   o get the total number of arrays
      iLenMOIO = _MAX_IO_GRPS * _MAX_IO_FAMS

c   o initialize the lists common block
      call getrec(-1,'JOBARC','MOIOVEC',iLenMOIO,moio)
      call getrec(-1,'JOBARC','MOIOWRD',iLenMOIO,moiowd)
      call getrec(-1,'JOBARC','MOIOSIZ',iLenMOIO,moiosz)
      call getrec(-1,'JOBARC','MOIODIS',iLenMOIO,moiods)
      call getrec(-1,'JOBARC','MOIOFIL',iLenMOIO,moiofl)
      call getrec(-1,'JOBARC','MOIOMXSZ',iLenMOIO,moiomxsz)

c   o initialize the distribution types
      call getrec(-1,'JOBARC','ISYMTYP',2*_MAX_IO_FAMS,isytyp)

c   o initialize the io_ptrs common block
      call getrec(-1,'JOBARC','TOTRECMO',_MAX_IO_LUNS,pRec)
      call getrec(-1,'JOBARC','TOTWRDMO',_MAX_IO_LUNS,iIntOff)
      if (pRec(1).eq.0) then
         do ihFile = 1, _MAX_IO_LUNS
            pRec(ihFile) = 1
         end do
      end if

c   o initialize the auxiliary cache
      call izero(quikget,iLenMOIO)

c   o initialize the automatic cache
      if (bInitCache) then

c        (sets /filspc/iprcln needed by aces_io_open)
         nprcwd = iflags(37)
         nslots = iflags(44)
         call aces_cache_init(iCore,i0,iMem,nprcwd,nslots)

#ifndef OPEN_LISTS_AS_NEEDED
c      o open the list files
         do ihFile = 1, _MAX_IO_LUNS
            iUnit = _FIRST_IO_LUN-1 + ihFile
            iPos = iieq(iLenMOIO,moiofl,1,iUnit)
            if (iPos.ne.0) then
               call gfname(szFiles(ihFile),szFile,iLength)
               inquire(file=szFile(1:iLength),exist=bExist,
     &                 err=666,iostat=iStat)
               if (bExist) then
                  call aces_io_open(ihFile)
               else
                  print *, '@ACES_IO_INIT: "',szFile(1:iLength),
     &                     '" - file not found'
                  call aces_exit(1)
               end if
            end if
         end do
#endif

c     end if (bInitCache)
      end if

      return

c   o I/O error
 666  print *, '@ACES_IO_INIT: I/O error'
      print '(/)'
      call aces_io_error('ACES_IO_INIT',iUnit,iStat)

c     end subroutine aces_io_init
      end

