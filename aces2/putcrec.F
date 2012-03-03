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

#define _MAXLEN 1024 /* 4096 or 8192 characters */
#define _IINTLN 8    /* PUTREC always write out 8-byte ints now. */ 

      subroutine putcrec(xFlag,szArchive,szRecName,iLength,szSrc)
      implicit none

c ARGUMENTS
      integer xFlag, iLength
      character*(*) szArchive, szRecName, szSrc(*)

c EXTERNAL FUNCTIONS
      integer iachar

c INTERNAL VARIABLES
      integer i, iNdx, nLeft, iRecLen
      integer*8 iBuf(_MAXLEN)

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert substring fits in iBuf
      if (iLength.gt.(_MAXLEN*_IINTLN)) then
         print *, '@PUTCREC: Assertion failed.'
         print *, '   iLength = ',iLength
         print *, '   maximum = ',_MAXLEN*_IINTLN
         i = 1
      end if
cYAU: Since szSrc is now an array pointer, we have no idea if iLength is good.
cc   o assert substring is within szSrc
c      if (len(szSrc).lt.iLength) then
c         print *, '@PUTCREC: Assertion failed.'
c         print *, '   iLength = ',iLength
c         print *, '   len(sz) = ',len(szSrc)
c         i = 1
c      end if
c   o assert ichar and iachar are the same
      if (ichar(' ').ne.iachar(' ')) then
         print *, '@PUTCREC: Assertion failed.'
         print *, '   iachar(szSpace) = ',iachar(' ')
         print *, '    ichar(szSpace) = ',ichar(' ')
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o iRecLen is the number of whole integers that contain szSrc(1:iLength)
      iRecLen = (iLength+_IINTLN-1)/_IINTLN

c   o flush the last integer with spaces
      if (iand(iLength,_IINTLN-1).ne.0) then
         i = ichar(' ')
         call c_memset(iBuf(iRecLen),i,_IINTLN)
      end if

c   o copy the string to the integer buffer
      call c_memmove(iBuf,szSrc,iLength)

c   o write to the job archive
      call putrec(xFlag,szArchive,szRecName,iRecLen,iBuf)

      return
c     end subroutine putcrec
      end

