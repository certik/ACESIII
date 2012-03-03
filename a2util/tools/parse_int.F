
c This subroutine scans a string and returns the first valid integer
c that falls into one of the following two regular expressions:
c         " *[+-]? *[0-9]+"     -->  1,   +1,   -1
c     " *\( *[+-]? *[0-9]+.*\)" --> (1), (+1), (-1)
c If the string is blank, then iVal=UNDEFINED, nChars=0, iErrPos=0.

c INPUT
c char*(*) sz : the string to parse

c OUTPUT
c int iVal    : the valid integer (undefined if iErrPos!=0)
c int nChars  : the number of characters processed (undefined if iErrPos!=0)
c int iErrPos : the position of an invalid character

#include "iachar.h"

      subroutine parse_int(sz,iVal,nChars,iErrPos)
      implicit none

c ARGUMENTS
      character*(*) sz
      integer iVal, nChars, iErrPos

c INTERNAL VARIABLES
      integer maxlen, i, j
      logical bParen, bNeg
      character*1 czSpace, czTab

c INTRINSIC FUNCTIONS
      integer   index
      intrinsic index
      character*1 achar
      intrinsic   achar

c ----------------------------------------------------------------------

      iErrPos = 0
      maxlen  = len(sz)
      czSpace = achar(_IACHAR_SPACE)
      czTab   = achar(_IACHAR_TAB)

      i = 1
      do while ((sz(i:i).eq.czSpace.or.sz(i:i).eq.czTab).and.
     &          i.le.maxlen)
         i = i + 1
      end do
      if (i.gt.maxlen) then
c      o the string is blank
         nChars = 0
         return
      end if

c   o parenthetical?
      bParen = (sz(i:i).eq.'(')
      if (bParen) then
         j = i + 1
         do while (sz(j:j).ne.')'.and.j.le.maxlen)
            j = j + 1
         end do
c      o no closing parenthesis
         if (j.gt.maxlen) iErrPos = i
         maxlen = j - 1
c      o find the next non-blank character
         i = i + 1
         do while ((sz(i:i).eq.czSpace.or.sz(i:i).eq.czTab).and.
     &             i.le.maxlen)
            i = i + 1
         end do
c      o the string is blank
         if (i.gt.maxlen) iErrPos = i
         if (iErrPos.ne.0) return
      end if

c   o sign?
      j = index('0123456789-+',sz(i:i))
      if (j.eq.0) then
         iErrPos = i
         return
      end if
      bNeg = (j.eq.11)
      if (j.gt.10) then
c      o find the next non-blank character
         j = i + 1
         do while ((sz(j:j).eq.czSpace.or.sz(j:j).eq.czTab).and.
     &             j.le.maxlen)
            j = j + 1
         end do
c      o the string is missing the number
         if (j.gt.maxlen) then
            iErrPos = i
            return
         end if
         i = j
      end if

c   o process all integers from i to maxlen
      iVal = 0
      j = index('0123456789',sz(i:i))
      if (j.eq.0) then
c      o no first number
         iErrPos = i
         return
      end if
      do while (j.ne.0.and.i.le.maxlen)
         iVal = 10*iVal + (j-1)
         i = i + 1
         j = index('0123456789',sz(i:i))
      end do
      if (bNeg) iVal = -iVal

c   o record the number of characters processed
      nChars = i - 1
      if (bParen) nChars = maxlen+1

c ----------------------------------------------------------------------

      return
c     end subroutine parse_int
      end

