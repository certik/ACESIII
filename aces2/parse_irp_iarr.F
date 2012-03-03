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

c This routine parses a string, szData, and extracts integers for nIrrps
c up to dim1 in groups of up to dim2. If nIrrps is zero, then it is set to
c the number of irreps in the first group up to dim1. The elements are
c separated by dashes, elements can be enclosed by parentheses, and
c groups are separated by forward slashes. nSpins is returned, telling
c the caller how many groups were read. iArr is initialized to zero.

c EXAMPLE:
c    szData=' 1 - 2 - 3 - 4 / (-5) - (-6) - 7 - 8 '
c This string yields: { {1,2,3,4} {-5,-6,7,8} }

      subroutine parse_irp_iarr(szData,nIrrps,nSpins,iArr,dim1,dim2)
      implicit none

c ARGUMENTS
      character*(*) szData
      integer nIrrps, nSpins, dim1, dim2
      integer iArr(dim1,dim2)

c EXTERNAL FUNCTIONS
      integer fnblnk

c INTERNAL VARIABLES
      integer iIrp, iSpin, i, j
      integer iLength, pStart, iFirst, nChars, iErrPos
      character*(1) czData
      logical bDone

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iIrp = 0
c   o assert nIrrps is in [0,inf)
      if (nIrrps.lt.0) then
         print *, '@PARSE_IRP_IARR: The expected number of irreps ',
     &            'cannot be negative.'
         print *, '   nIrrps = ',nIrrps
         iIrp = -1
      end if
c   o assert dim1 and dim2 are natural
      if ((dim1.lt.1).or.(dim2.lt.1)) then
         print *, '@PARSE_IRP_IARR: The destination array is ',
     &            'ill-defined.'
         print *, '   dim1 = ',dim1
         print *, '   dim2 = ',dim2
         iIrp = -1
      end if
c   o assert dim1 is not less than nIrrps
      if (dim1.lt.nIrrps) then
         print *, '@PARSE_IRP_IARR: The destination array is ',
     &            'too small for the expected data.'
         print *, '   nIrrps = ',nIrrps
         print *, '   dim1   = ',dim1
         iIrp = -1
      end if
      if (iIrp.ne.0) then
         call errex
      end if
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o initialize the destination array
      do j = 1, dim2
      do i = 1, dim1
         iArr(i,j) = 0
      enddo
      enddo
      

c   o quit if the string is empty
      iFirst = fnblnk(szData)
      if (iFirst.eq.0) then
         nSpins = 0
         return
      end if

c   o initialize the iArr indices
      iIrp  = 1
      iSpin = 1

c   o start reading at the beginning and continue until nothing is left
      pStart  = iFirst
      iLength = len(szData)
      bDone   = .false.
      do while (.not.bDone)

c      o find a number
         nChars = 0
         if (pStart.le.iLength) then
            call parse_int(szData(pStart:),iArr(iIrp,iSpin),
     &                     nChars,iErrPos)
            if (iErrPos.ne.0) then
               print *, '@PARSE_SET_IARR: error near position ',
     &                  pStart - 1 + iErrPos
               print *, '"',szData,'"'
               call errex
            end if
            pStart = pStart + nChars
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

         if (iFirst.ne.0) then
c         o process more text

            i = pStart - 1 + iFirst
            czData = szData(i:i)
            if (czData.ne.'-'.and.czData.ne.'/') then
               print *, '@PARSE_IRP_IARR: invalid delimiter at ',
     &                  'position ',i
               print *, '"',szData,'"'
               call errex
            end if
            pStart = i + 1

            if (czData.eq.'-') then
c            o read another irrep value
               if (iIrp.eq.dim1) then
                  print *, '@PARSE_IRP_IARR: number of irreps exceeds ',
     &                     'buffer near position ',pStart
                  print *, '"',szData,'"'
                  call errex
               end if
               iIrp = iIrp + 1
            else
c            o start a new spin column
               if (iSpin.eq.dim2) then
                  print *, '@PARSE_IRP_IARR: too many spin sets'
                  print *, '"',szData,'"'
                  call errex
               end if
               iSpin = iSpin + 1
               if (nIrrps.eq.0) nIrrps = iIrp
               if (iIrp.ne.nIrrps) then
                  print *, '@PARSE_IRP_IARR: number of irreps is ',
     &                     'inconsistent.'
                  print *, '                 ',nIrrps,' expected',
     &                     ' but read ',iIrp,' near position ',pStart
                  print *, '"',szData,'"'
                  call errex
               end if
               iIrp = 1
            end if

c        else if (iFirst.eq.0) then
         else
c         o no text left

            if (nIrrps.eq.0) nIrrps = iIrp
            if (iIrp.ne.nIrrps) then
               print *, '@PARSE_IRP_IARR: number of irreps is ',
     &                  'inconsistent.'
               print *, '                 ',nIrrps,' expected',
     &                  ' but read ',iIrp,' near position ',pStart
               print *, '"',szData,'"'
               call errex
            end if
            bDone = .true.

c        end if (iFirst)
         end if

c     end do while (.not.bDone)
      end do

c   o record the output
      nSpins = iSpin

      return
      end

