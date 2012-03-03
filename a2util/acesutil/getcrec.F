
c This routine reads a character string from the job archive.

c WARNING!
c    The "record length" that is associated with the logical record is the
c number of whole integers that contain the string. Therefore, calling getrec
c with IFLAG=0 may not return the actual length of the string.

c INPUT
c int      IFLAG     : (same as getrec)
c char*(*) SZARCHIVE : (same as getrec)
c char*(*) SZRECNAME : (same as getrec)

c INPUT/OUTPUT
c int ILENGTH : on input, this is the substring-length of the record to get
c               on output, if IFLAG=0, the string-length of the logical record
c                  is returned

c OUTPUT
c char*(*) SZDEST : (same as getrec)

#define _MAXLEN 1024 /* 4096 or 8192 characters */

#ifdef F_64BIT
#  define _IINTLN 8
#else
#  define _IINTLN 4
#endif

      subroutine getcrec(iFlag,szArchive,szRecName,iLength,szDest)
      implicit none

c ARGUMENTS
      integer iFlag, iLength
      character*(*) szArchive, szRecName, szDest(*)

c INTERNAL VARIABLES
      integer i, iNdx, nLeft, iBuf(_MAXLEN), iRecLen

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert substring fits in iBuf
      if (iLength.gt.(_MAXLEN*_IINTLN)) then
         print *, '@GETCREC: Assertion failed.'
         print *, '   iLength = ',iLength
         print *, '   maximum = ',_MAXLEN*_IINTLN
         i = 1
      end if
cYAU: Since szDest is now an array pointer, we have no idea if iLength is good.
cc   o assert substring is within szDest
c      if (len(szDest).lt.iLength) then
c         print *, '@GETCREC: Assertion failed.'
c         print *, '   iLength = ',iLength
c         print *, '   len(sz) = ',len(szDest)
c         i = 1
c      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o assume the length of the record from iLength
      iRecLen = (iLength+_IINTLN-1)/_IINTLN

c   o read from the job archive
      call getrec(iFlag,szArchive,szRecName,iRecLen,iBuf)

c   o comply with getrec queries
      if (iFlag.eq.0) then
         iLength = iRecLen*_IINTLN
         return
      end if

c   o copy the string from the integer buffer
      call c_memmove(szDest,iBuf,iLength)

      return
c     end subroutine getcrec
      end

