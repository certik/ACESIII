
c This routine removes all traces of a storage unit.

c INPUT
c int      IUNIT  : the external file unit to destroy
c char*(*) SZFILE : the internal name of the file to destroy
c                   NOTE: This is only used if the external unit is closed.

c#define _DEBUG_ACES_IO_REMOVE

      subroutine aces_io_remove(iUnit,szFile)
      implicit none

c ARGUMENTS
      integer iUnit
      character*(*) szFile

c INTERNAL VARIABLES
      integer i, iStat
      character*80 szExt
      integer      iLength
      integer ihFile, iFamNdx, iStorLen
      logical bExist, bOpened

c COMMON BLOCKS
#include "lists.com"
#include "cache.com"
#include "sympop.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_IO_REMOVE: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         i = 1
      end if
c   o assert iUnit is properly bound
      if ((iUnit.lt._FIRST_IO_LUN).or.
     &    (_FIRST_IO_LUN-1+_MAX_IO_LUNS.lt.iUnit)) then
         print *, '@ACES_IO_REMOVE: Assertion failed.'
         print *, '   iUnit = ',iUnit
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o delete the unit and/or file
      inquire(unit=iUnit,opened=bOpened,err=666,iostat=iStat)
      if (bOpened) then
#ifdef _DEBUG_ACES_IO_REMOVE
         print *, '@ACES_IO_REMOVE: deleting unit ',iUnit,' / ',szFile
#endif
         close(unit=iUnit,status='DELETE',err=666,iostat=iStat)
c      o destroy the relevant cache slots
         do i = 1, cachnum
            if (cachfil(i).eq.iUnit) then
               cachfil(i)  = 0
               cachrec(i)  = 0
               cachmod(i)  = 0
               lrustats(i) = 0
            end if
         end do
      else
         call gfname(szFile,szExt,iLength)
         inquire(file=szExt(1:iLength),exist=bExist,
     &           err=666,iostat=iStat)
         if (bExist) then
#ifdef _DEBUG_ACES_IO_REMOVE
            print *, '@ACES_IO_REMOVE: deleting file ',szExt(1:iLength)
#endif
#ifdef _UNICOS
            call pxfunlink(szExt(1:iLength),iLength,iStat)
            if (iStat.ne.0) goto 666
#else
            call f_remove(szExt(1:iLength))
#endif
         end if
      end if

c   o get the internal unit index and the first family index
      ihFile  = 1 + iUnit - _FIRST_IO_LUN
      iFamNdx = 1 + (ihFile-1)*_FAMS_PER_FILE
#ifdef _DEBUG_ACES_IO_REMOVE
      print *, '@ACES_IO_REMOVE: deleting list families ',
     &         iFamNdx,' through ',iFamNdx-1+_FAMS_PER_FILE
#endif

c   o reset the list stats
      iStorLen = _MAX_IO_GRPS*_FAMS_PER_FILE
      call izero(moio  (1,iFamNdx),iStorLen)
      call izero(moiowd(1,iFamNdx),iStorLen)
      call izero(moiods(1,iFamNdx),iStorLen)
      call izero(moiosz(1,iFamNdx),iStorLen)
      call izero(moiofl(1,iFamNdx),iStorLen)
      call izero(moiomxsz(1,iFamNdx),iStorLen)

cYAU - There is so much bad programming related to this array that we cannot
c      destroy it.
cc   o reset the distribution types
c      call izero(isytyp(1,iFamNdx),2*_FAMS_PER_FILE)

c   o reset the free-space pointers
      pRec(ihFile)    = 1
      iIntOff(ihFile) = 0

c   o mark MOIO as modified
      bIOMod = .true.

      return

c   o I/O error
 666  print *, '@ACES_IO_REMOVE: I/O error'
      print *, '                 file = "',szFile,'"'
      print '(/)'
      call aces_io_error('ACES_IO_REMOVE',iUnit,iStat)

c     end subroutine aces_io_remove
      end

