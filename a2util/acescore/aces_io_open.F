
c This routine opens the unformatted, direct access files used by ACES
c for double precision array (list) storage.

c INPUT
c int IHFILE : the internal unit number (handle) of the file to open

c#define _DEBUG_ACES_IO_OPEN

      subroutine aces_io_open(ihFile)
      implicit none

c ARGUMENTS
      integer ihFile

c INTERNAL VARIABLES
#include "szfiles.dat"
      character*80 szFile
      integer      iLength
      integer iUnit, iStat
      logical bExist

c COMMON BLOCKS
#include "lists.com" /* for bIOUp */
#include "filspc.com" /* for iprcrl */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iStat = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_IO_OPEN: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iStat = 1
      end if
c   o assert iprcrl was initialized
      if (iprcrl.lt.1) then
         print *, '@ACES_IO_OPEN: Assertion failed.'
         print *, '   iprcrl = ',iprcrl
         iStat = 1
      end if
c   o assert ihFile is properly bound
      if ((ihFile.lt.1).or.(_MAX_IO_LUNS.lt.ihFile)) then
         print *, '@ACES_IO_OPEN: Assertion failed.'
         print *, '   ihFile = ',ihFile
         iStat = 1
      end if
      if (iStat.ne.0) call aces_exit(iStat)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      call gfname(szFiles(ihFile),szFile,iLength)
      inquire(file=szFile(1:iLength),exist=bExist,err=666,iostat=iStat)

#ifdef _DEBUG_ACES_IO_OPEN
      print *, '@ACES_IO_OPEN: unit=',ihFile-1+_FIRST_IO_LUN,
     &                      ', file=',szFile(1:iLength),
     &                      ', recl=',iprcrl
      print *, '               (iprcwd=',iprcwd,')'
#endif

      iUnit = ihFile-1+_FIRST_IO_LUN
      if (bExist) then
         open(unit=iUnit,file=szFile(1:iLength),
     &        form='UNFORMATTED',access='DIRECT',recl=iprcrl,
     &        status='OLD',err=666,iostat=iStat)
      else
         open(unit=iUnit,file=szFile(1:iLength),
     &        form='UNFORMATTED',access='DIRECT',recl=iprcrl,
     &        status='NEW',err=666,iostat=iStat)
         write(unit=iUnit,rec=1,err=666,iostat=iStat) 0
      end if
      return

 666  call aces_io_error('ACES_IO_OPEN',ihFile-1+_FIRST_IO_LUN,iStat)

c     end subroutine aces_io_open
      end

