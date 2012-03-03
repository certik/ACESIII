
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
      integer      iLen
      integer iStat
      logical bExist

c COMMON BLOCKS
#include "lists.com" /* for bIOUp */
#include "filspc.com" /* for iprcrl */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iLen = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_IO_OPEN: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iLen = 1
      end if
c   o assert iprcrl was initialized
      if (iprcrl.lt.1) then
         print *, '@ACES_IO_OPEN: Assertion failed.'
         print *, '   iprcrl = ',iprcrl
         iLen = 1
      end if
c   o assert ihFile is properly bound
      if ((ihFile.lt.1).or.(_MAX_IO_LUNS.lt.ihFile)) then
         print *, '@ACES_IO_OPEN: Assertion failed.'
         print *, '   ihFile = ',ihFile
         iLen = 1
      end if
      if (iLen.ne.0) call aces_exit(iLen)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      call gfname(szFiles(ihFile),szFile,iLen)
      inquire(file=szFile(1:iLen),exist=bExist)

#ifdef _DEBUG_ACES_IO_OPEN
      print *, '@ACES_IO_OPEN: unit=',ihFile-1+_FIRST_IO_LUN,
     &                      ', file=',szFile(1:iLen),
     &                      ', recl=',iprcrl
      print *, '               (iprcwd=',iprcwd,')'
#endif

      if (bExist) then
         open(unit=ihFile-1+_FIRST_IO_LUN,file=szFile(1:iLen),
     &        form='UNFORMATTED',access='DIRECT',recl=iprcrl,
     &        status='OLD',err=9999,iostat=iStat)
      else
         open(unit=ihFile-1+_FIRST_IO_LUN,file=szFile(1:iLen),
     &        form='UNFORMATTED',access='DIRECT',recl=iprcrl,
     &        status='NEW',err=9999,iostat=iStat)
      end if
      return

 9999 call aces_io_error('ACES_IO_OPEN',ihFile-1+_FIRST_IO_LUN,iStat)

c     end subroutine aces_io_open
      end

