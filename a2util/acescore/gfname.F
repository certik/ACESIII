
c This routine maps an internal file name (limited to 8 characters)
c like 'JOBARC' to a user-defined (external) filename like '../JOBARC.old'.

cjp - For MR-BWCC, various files are affixed with ".##", in which the hashes
c     represent the 2-digit process id.
cYAU - Our parallel file processing has not been fully designed. Once it is, I
c      imagine this 'feature' will be refined or replaced.

c INPUT
c char*(*) SZINT : the internal file name

c OUTPUT
c char*80 SZEXT   : the external file name for use with OPEN or INQUIRE
c int     ILENGTH : the character length of SZEXT

c#define _DEBUG_GFNAME
#define _FILES_UNIT 92 /* the file unit number for FILES (temporary) */

      subroutine gfname(szInt,szExt,iLength)
      implicit none

c ARGUMENTS
      character*(*) szInt, szExt
      integer iLength

c EXTERNAL FUNCTIONS
      integer fnblnk
      character*1 achar

c PARAMETERS
      character*1 czSpace
      parameter (czSpace=' ')

c INTERNAL VARIABLES
      character*80 szTmp
      integer iEnd, iEnd1, iTmp, iStat
      logical bExist, bDone

c COMMON BLOCKS
#include "bwcc.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert szInt is not empty
      if (fnblnk(szInt).eq.0) then
         print *, '@GFNAME: Assertion failed.'
         print *, '   szInt = "',szInt,'"'
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o determine the shortest string length of szInt
      iTmp = min(8,len(szInt))
      iEnd = 0
      do while ((szInt(iEnd+1:iEnd+1).ne.czSpace).and.(iEnd.lt.iTmp))
         iEnd = iEnd + 1
      end do

c   o initialize szExt
      if (iEnd.gt.0) then
         szExt(1:iEnd) = szInt(1:iEnd)
         iLength = iEnd
      else
         iLength = 0
         return
      end if

c   o attempt to map the filename to a user-defined file
      inquire(file='FILES',exist=bExist,err=666,iostat=iStat)
      if (bExist) then
         open(unit=_FILES_UNIT,file='FILES',
     &        status='OLD',form='FORMATTED',
     &        err=666,iostat=iStat)
         rewind(_FILES_UNIT,err=666,iostat=iStat)
         bDone = .false.
         do while (.not.bDone)
            read(unit=_FILES_UNIT,fmt='(a)',
     &           end=200,err=666,iostat=iStat) szTmp
            if ((szTmp(1:iEnd).eq.szInt(1:iEnd)) .and.
     &          (szTmp(iEnd+1:iEnd+1).eq.czSpace)     ) then
               iEnd1 = iEnd+1
               do while ((szTmp(iEnd1+1:iEnd1+1).ne.czSpace).and.
     &                   (iEnd1.lt.80))
                  iEnd1 = iEnd1 + 1
               end do
               iLength = iEnd1-(iEnd+1)
               szExt(1:iLength) = szTmp(iEnd+2:iEnd1)
               bDone = .true.
            end if
         end do
 200     continue
         close(_FILES_UNIT,status='KEEP',err=666,iostat=iStat)
c     end if (bExist)
      end if

c   o add the process id to the end of various files
      if (.not.masik.and.isbwcc) then
         if (szInt(1:iEnd).eq.'MOINTS'.or.
     &       szInt(1:iEnd).eq.'MOABCD'.or.
     &       szInt(1:iEnd).eq.'FOCKCD'    ) then
            szExt(iLength+1:iLength+1) = '.'
            szExt(iLength+2:iLength+2) = achar(48+iref/10)
            szExt(iLength+3:iLength+3) = achar(48+iref-10*(iref/10))
            iLength = iLength+3
         end if
      end if

#ifdef _DEBUG_GFNAME
      print *, '@GFNAME: (debug) "',szInt,'" -> "',szExt(1:iLength),'"'
#endif

      return

c   o FILES I/O error
 666  print *, '@GFNAME: I/O error on FILES'
      print *, '         internal file name = "',szInt,'"'
      print '(/)'
      call aces_io_error('GFNAME',_FILES_UNIT,iStat)

c     end subroutine gfname
      end

