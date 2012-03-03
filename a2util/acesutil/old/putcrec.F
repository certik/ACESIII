
c This routine writes a character string to the job archive.

c WARNING!
c    The "record length" that gets associated with the logical record is the
c number of whole integers that contain the string. Therefore, calling getrec
c with the query flag may not return the actual length of the string.

c INPUT
c int      XFLAG     : (same as putrec)
c char*(*) SZARCHIVE : (same as putrec)
c char*(*) SZRECNAME : (same as putrec)
c int      ILENGTH   : the substring-length of the record to store
c char*(*) SZSRC     : (same as putrec)

#define _MAXLEN 256 /* 1024 or 2048 characters */

#ifdef F_64BIT
#  define _IINTLN 8
#else
#  define _IINTLN 4
#endif

      subroutine putcrec(xFlag,szArchive,szRecName,iLength,szSrc)
      implicit none

c ARGUMENTS
      integer xFlag, iLength
      character*(*) szArchive, szRecName, szSrc

c INTERNAL VARIABLES
      integer i, iNdx, nLeft, iBuf(_MAXLEN), iRecLen

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert substring fits in iBuf
      if (iLength.gt.(_MAXLEN*_IINTLN)) then
         print *, '@PUTCREC: Assertion failed.'
         print *, '   iLength = ',iLength
         print *, '   maximum = ',_MAXLEN*_IINTLN
         i = 1
      endif
c   o assert substring is within szSrc
      if (len(szSrc).lt.iLength) then
         print *, '@PUTCREC: Assertion failed.'
         print *, '   iLength = ',iLength
         print *, '   len(sz) = ',len(szSrc)
         i = 1
      endif
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o iNdx is the moving szSrc index
      iNdx = 1

c   o iRecLen is the number of whole integers that contain szSrc(1:iLength)
c     (residual handled later)
      iRecLen = iLength/_IINTLN

c   o the format definition for reading whole integers from szSrc
 10   format(a _IINTLN)

c   o read from szSrc into iBuf by whole integers
      if (iRecLen.gt.0) then
         do i = 1, iRecLen
            read(szSrc(iNdx:iNdx+_IINTLN-1),10) iBuf(i)
            iNdx = iNdx + _IINTLN
         end do
      end if

c   o pick up the slack (while flushing the integer with spaces)
      nLeft = iand(iLength,_IINTLN-1)
      if (nLeft.ne.0) then
         iRecLen = iRecLen + 1
         read(szSrc(iNdx:iNdx+nLeft-1),10) iBuf(iRecLen)
      end if

c   o write to the job archive
      call putrec(xFlag,szArchive,szRecName,iRecLen,iBuf)

      return
c     end subroutine putcrec
      end

