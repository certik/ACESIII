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

c This routine parses a string, szData, and extracts a one-dimensional
c array of integers up to nMax. The elements are separated by forward slashes,
c ranges are defined with dashes or '>', and elements can be enclosed by
c parentheses. nInts is returned, telling the caller how many integers were
c read. iArr is NOT initialized to zero.

c EXAMPLE:
c    szData=' 1 / 3 / 4>6 / (-7)-(-9) '
c This string yields: { 1, 3, 4, 5, 6, -7, -8, -9 }

c NOTES:
c  - The data is not checked for duplicates; therefore, '1/2/1-3' will
c    yield 5 separate elements.
c  - The data is not checked for consecutiveness, meaning the caller
c    must assume the data is in random order.

      subroutine parse_set_iarr(szData,nInts,iArr,nMax)
      implicit none

c ARGUMENTS
      character*(*) szData
      integer nInts, nMax
      integer iArr(nMax)

c EXTERNAL FUNCTIONS
      integer fnblnk

c INTERNAL VARIABLES
      integer iInts
      integer iStart, iEnd, inc, nChars, iErrPos
      integer iLength, pStart, iFirst, i
      character*(1) czData
      logical bDone

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iInts = 0
c   o assert nMax is natural
      if (nMax.lt.1) then
         print *, '@PARSE_SET_IARR: The destination array is ',
     &            'ill-defined.'
         print *, '   nMax = ',nMax
         iInts = -1
      end if
      if (iInts.ne.0) then
         call errex
      end if
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o quit if the string is empty
      iFirst = fnblnk(szData)
      if (iFirst.eq.0) then
         nInts = 0
         return
      end if

c   o initialize the iArr index (we increment just before writing a
c     new element)
      iInts = 0

c   o start reading from the first character
      pStart  = iFirst
      iLength = len(szData)
      bDone   = .false.
      do while (.not.bDone)

c      o find a number
         nChars = 0
         if (pStart.le.iLength) then
            call parse_int(szData(pStart:),iStart,nChars,iErrPos)
            if (iErrPos.ne.0) then
               print *, '@PARSE_SET_IARR: error near position ',
     &                  pStart - 1 + iErrPos
               print *, '"',szData,'"'
               call errex
            end if
            pStart = pStart + nChars
            iEnd   = iStart
         end if
         if (nChars.eq.0) then
            print *, '@PARSE_SET_IARR: missing value at position ',
     &               pStart
            print *, '"',szData,'"'
            call errex
         end if

c      o point to the next character
         iFirst = 0
         if (pStart.le.iLength) iFirst = fnblnk(szData(pStart:))

c      o check for a range
         if (iFirst.ne.0) then
            i = pStart - 1 + iFirst
            czData = szData(i:i)
            if (czData.eq.'-'.or.czData.eq.'>') then
               nChars = 0
               if (i.lt.iLength) then
                  call parse_int(szData(i+1:),iEnd,nChars,iErrPos)
                  if (iErrPos.ne.0) then
                     print *, '@PARSE_SET_IARR: error near position ',
     &                        i + iErrPos
                     print *, '"',szData,'"'
                     call errex
                  end if
                  pStart = i + nChars + 1
               end if
               if (nChars.eq.0) then
                  print *, '@PARSE_SET_IARR: missing end of range ',
     &                     'after position ',i
                  print *, '"',szData,'"'
                  call errex
               end if
            end if
         end if

         inc = 1
         if (iStart.gt.iEnd) inc = -1

c      o check for overflow
         i = abs(iEnd+inc-iStart)
         if (iInts+i.gt.nMax) then
            print *, '@PARSE_SET_IARR: exceeded limit of',nMax,' values'
            print *, '"',szData,'"'
            call errex
         end if

c      o enter the values
         do i = iStart, iEnd, inc
            iInts = iInts + 1
            iArr(iInts) = i
         end do

c      o point to the next character
         iFirst = 0
         if (pStart.le.iLength) iFirst = fnblnk(szData(pStart:))
         bDone = (iFirst.eq.0)

         if (.not.bDone) then
            i = pStart - 1 + iFirst
            czData = szData(i:i)
            if (czData.ne.'/') then
               print *, '@PARSE_SET_IARR: invalid character at ',
     &                  'position ',i
               print *, '"',szData,'"'
               call errex
            end if
            pStart = i + 1
         end if

c     end do while (.not.bDone)
      end do

c      print *, 'DEBUG: ',(iArr(i),i=1,iInts)

c   o record the output
      nInts = iInts

      return
      end

