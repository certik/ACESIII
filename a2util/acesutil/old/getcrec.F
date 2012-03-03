
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

#define _MAXLEN 256 /* 1024 or 2048 characters */

#ifdef F_64BIT
#  define _IINTLN 8
#else
#  define _IINTLN 4
#endif

      subroutine getcrec(iFlag,szArchive,szRecName,iLength,szDest)
      implicit none

c ARGUMENTS
      integer iFlag, iLength
      character*(*) szArchive, szRecName, szDest

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
      endif
c   o assert substring is within szDest
      if (len(szDest).lt.iLength) then
         print *, '@GETCREC: Assertion failed.'
         print *, '   iLength = ',iLength
         print *, '   len(sz) = ',len(szDest)
         i = 1
      endif
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

c   o iNdx is the moving szDest index
      iNdx = 1

c   o iRecLen is the number of whole integers that contain szDest(1:iLength)
c     (residual handled later)
      iRecLen = iLength/_IINTLN

c   o the format definition for reading whole integers into szDest
 10   format(a _IINTLN)

c   o read from iBuf into szDest by whole integers
      if (iRecLen.gt.0) then
         do i = 1, iRecLen
            write(szDest(iNdx:iNdx+_IINTLN-1),10) iBuf(i)
            iNdx = iNdx + _IINTLN
         end do
      end if

c   o pick up the slack
      nLeft = iand(iLength,_IINTLN-1)
      if (nLeft.ne.0) then
         iRecLen = iRecLen + 1
         do i = 1, nLeft
            write(szDest(iNdx:iNdx),'(a1)') iBuf(iRecLen)
            iBuf(iRecLen) = rshift(iBuf(iRecLen),8)
            iNdx = iNdx + 1
         end do
      end if

      return
c     end subroutine getcrec
      end

