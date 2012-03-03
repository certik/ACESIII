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
c
c      program main
c      implicit none
c      integer iErr
c      open(unit=13,file='CHSSI.MOL',form='formatted',status='unknown')
c      close(13,status='delete')
c      call segment_mol('MOL','CHSSI.MOL',12,13,iErr)
c      if (iErr.ne.0) stop 'Failed to segment MOL file.'
c      print '(a)', 'Successfully created CHSSI.MOL'
c      end

c ----------------------------------------------------------------------

c This routine reads a MOL file and segments the basis into (mostly)
c disjoint subshells.

c INPUT:
c char*(*) szOldMOL : file name in MOL format to convert
c char*(*) szNewMOL : file name in MOL format to create
c int      iUOld    : unit number to use for old file
c int      iUNew    : unit number to use for new file
c int      iErr     : return code (success==0, failure!=0)

      subroutine segment_mol(szOldMOL,szNewMOL,iUOld,iUNew,iErr)
      implicit none

c   o ARGUMENTS
      character*(*) szOldMOL, szNewMOL
      integer iUOld, iUNew, iErr

c   o PARAMETERS
      integer maxLines, maxAngMom
      parameter (maxLines=10000, maxAngMom=8)
      integer maxPpB, maxCpB
      parameter (maxPpB=500, maxCpB=500)

c   o VARIABLES
      double precision dTmp, dCharge
      double precision dPrim(maxPpB), dCoef(maxPpB,maxCpB)
      integer iAtom, iLine, iOldLine, iNewLine
      integer iAngMom, nAMShells
      integer nOldBpAM(1+maxAngMom), iOldBlock, nOldBlocks, nBlocksLeft
      integer nNewBpAM(1+maxAngMom), iNewBlock, nNewBlocks
      integer iPrim, nOldPrims, nNewPrims
      integer iCont, nOldConts, nNewConts
      integer iStart(maxCpB), iEnd(maxCpB), iStartPrv, iEndPrv
      integer nSubShells, SubShellPtr(maxCpB)
      integer i, j, iOff, iReg0, iReg1, iReg2, iReg3
      character*80 szOldBuf(maxLines), szNewBuf(maxLines)
      character*80 szError, szLine
      logical bExist

c ----------------------------------------------------------------------

c   o open the files or return
      iErr = 0
c      inquire(file=szNewMOL,exist=bExist)
c      if (bExist) return
      szError = '@SEGMENT_MOL: Error opening '//szOldMOL
      open(unit=iUOld,file=szOldMOL,form='formatted',status='old',
     &     err=666,iostat=iErr)
      szError = '@SEGMENT_MOL: Error opening '//szNewMOL
      open(unit=iUNew,file=szNewMOL,form='formatted',
     &     err=666,iostat=iErr)

c ----------------------------------------------------------------------

c   o create the 5-line header
      szError = '@SEGMENT_MOL: Error copying header'
      do i = 1, 5
         read( iUOld,'(a)',err=666,iostat=iErr) szLine
         write(iUNew,'(a)',err=666,iostat=iErr) szLine
      end do

c   o loop over (symmetry unique) atoms
      szError = '@SEGMENT_MOL: Error reading new atom data'
      read(iUOld,'(a)',err=666,iostat=iErr) szOldBuf(1)
      iAtom = 0
      do while (szOldBuf(1)(1:6).ne.'FINISH')
         iAtom = iAtom + 1

c      o flush the file buffers
         iNewLine = 0
         szNewBuf(1) = ' '
         do iLine = 2, maxLines
            szOldBuf(iLine) = ' '
            szNewBuf(iLine) = ' '
         end do

c ----------------------------------------------------------------------

c   o read the block statistics from line 1
      read(szOldBuf(1),'(6x,f14.8,i5,10i5)')
     &    dCharge, iReg0, nAMShells
      if (nAMShells.gt.1+maxAngMom) then
         szError = '@SEGMENT_MOL: Angular momentum limit exceeded.'
         iErr = 1
         goto 666
      end if
      read(szOldBuf(1),'(6x,f14.8,i5,10i5)')
     &    dTmp, iReg0, iReg1, (nOldBpAM(i),i=1,nAMShells)

c   o count the old blocks and reset the new blocks
      nOldBlocks = 0
      do i = 1, nAMShells
         nOldBlocks  = nOldBlocks + nOldBpAM(i)
         nNewBpAM(i) = 0
      end do

c   o read the symbol, orbit index, and XYZ coordinates
      szError = '@SEGMENT_MOL: Error reading XYZ coordinates'
      read(iUOld,'(a)',err=666,iostat=iErr) szOldBuf(2)

c   o read in all of the shell blocks
      iOldBlock = 0
      iOldLine  = 2
      do while (iOldBlock.lt.nOldBlocks)
         iOldBlock = iOldBlock + 1

c      o read the number of primitives and contractions
         szError = '@SEGMENT_MOL: Error reading block header'
         iOldLine = iOldLine + 1
         read(iUOld,'(a)',err=666,iostat=iErr) szOldBuf(iOldLine)

c----------------------------------------------------------------------------
c   Change any null chars to blanks.
c----------------------------------------------------------------------------

         do i = 1, len(szOldBuf(iOldLine))
            if (szOldBuf(iOldLine)(i:i) .eq. char(0)) 
     *           szOldBuf(iOldLine)(i:i) = ' '
         enddo

         read(szOldBuf(iOldLine),'(2i5)') nOldPrims, nOldConts
         if (nOldPrims.gt.maxPpB.or.nOldConts.gt.maxCpB) then
            szError = '@SEGMENT_MOL: number of functions exceeded'
            iErr = 1
            goto 666
         end if

c      o load this block
         szError = '@SEGMENT_MOL: Error reading block data'
         do iPrim = 1, nOldPrims
            iOldLine = iOldLine + 1
            read(iUOld,'(a)',err=666,iostat=iErr) szOldBuf(iOldLine)
            do iCont = 4, nOldConts, 4
               iOldLine = iOldLine + 1
               read(iUOld,'(a)',err=666,iostat=iErr) szOldBuf(iOldLine)
            end do
         end do

c     end do (read old blocks)
      end do

c ----------------------------------------------------------------------

c   o write new block data into szNewBuf
      iAngMom = 0
      nBlocksLeft = nOldBpAM(1)
      iOldBlock = 0
      iOldLine  = 2
      do while (iOldBlock.lt.nOldBlocks)
         iOldBlock = iOldBlock + 1

c      o reset the primitive start and end indices for each contraction
         do i = 1, nOldConts
            iStart(i) = 0
            iEnd(i)   = 0
         end do

         iOldLine = iOldLine + 1
         read(szOldBuf(iOldLine),'(2i5)') nOldPrims, nOldConts
         do iPrim = 1, nOldPrims

c         o read the vector of exponents and array of coefficients
            iOldLine = iOldLine + 1
            read(szOldBuf(iOldLine),'(4f18.9)')
     &          dPrim(iPrim),(dCoef(iPrim,i),i=1,min(3,nOldConts))
            iCont = 4
            do while (iCont.le.nOldConts)
               iOldLine = iOldLine + 1
               read(szOldBuf(iOldLine),'(4f18.9)')
     &             (dCoef(iPrim,i),i=iCont,min(iCont+3,nOldConts))
               iCont = iCont + 4
            end do

c         o scan the coefficients and mark the primitives that each uses
            do iCont = 1, nOldConts
               if (dCoef(iPrim,iCont).ne.0.d0) then
                  if (iStart(iCont).eq.0) iStart(iCont) = iPrim
                  iEnd(iCont) = iPrim
               end if
            end do

         end do

c      o determine the number of new blocks coming from this old block
         nNewBlocks = 0
         iStartPrv  = 0
         iEndPrv    = 0
         do iCont = 1, nOldConts
            if (iStart(iCont).ne.iStartPrv.and.
     &          iEnd(iCont).ne.iEndPrv) then
               nNewBlocks = nNewBlocks + 1
               SubShellPtr(nNewBlocks) = iCont
               iStartPrv = iStart(iCont)
               iEndPrv   = iEnd(iCont)
            end if
         end do

c      o write the set of new blocks to the new file buffer
         do iNewBlock = 1, nNewBlocks
            nNewPrims = 1 + iEnd(SubShellPtr(iNewBlock))
     &                    - iStart(SubShellPtr(iNewBlock))
            if (iNewBlock.lt.nNewBlocks) then
               nNewConts =   SubShellPtr(iNewBlock+1)
     &                     - SubShellPtr(iNewBlock)
            else
               nNewConts = 1+nOldConts-SubShellPtr(iNewBlock)
            end if
c         o block header
            iNewLine = iNewLine + 1
            write(szNewBuf(iNewLine),'(2i5)') nNewPrims, nNewConts
            do iPrim = iStart(SubShellPtr(iNewBlock)),
     &                 iEnd(SubShellPtr(iNewBlock))
c            o exponent line
               iNewLine = iNewLine + 1
               write(szNewBuf(iNewLine),'(4f18.9)')
     &              dPrim(iPrim),
     &              (dCoef(iPrim,SubShellPtr(iNewBlock)+i),
     &               i=0,min(2,nNewConts-1))
c            o additional coefficient lines
               iOff = 3
               do while (iOff.lt.nNewConts)
                  iNewLine = iNewLine + 1
                  write(szNewBuf(iNewLine),'(4f18.9)')
     &                 (dCoef(iPrim,SubShellPtr(iNewBlock)+iOff+i),
     &                  i=0,min(3,nNewConts-1-iOff))
                  iOff = iOff + 4
               end do
            end do
c        end do iNewBlock = 1, nNewBlocks
         end do

c      o increment the angular momentum
         nNewBpAM(1+iAngMom) = nNewBpAM(1+iAngMom) + nNewBlocks
         nBlocksLeft = nBlocksLeft - 1
         if (nBlocksLeft.eq.0.and.iAngMom.lt.nAMShells) then
            iAngMom = iAngMom + 1
            nBlocksLeft = nOldBpAM(1+iAngMom)
         end if

c     end do (process old blocks)
      end do
c ----------------------------------------------------------------------

c      o write the new file buffer to the new file
         read(szOldBuf(1),'(6x,f14.8,i5,10i5)') dCharge, iReg0, iReg1
         szError = '@SEGMENT_MOL: Error writing new atom data'
         write(iUNew,'(6x,f14.8,i5,10i5)',err=666,iostat=iErr)
     &        dCharge, iReg0, iReg1, (nNewBpAM(i),i=1,nAMShells)
         write(iUNew,'(a)',err=666,iostat=iErr) szOldBuf(2)
         szError = '@SEGMENT_MOL: Error writing new basis set data'
         do iLine = 1, iNewLine
            write(iUNew,'(a)',err=666,iostat=iErr) szNewBuf(iLine)
         end do

c      o read the next basis set statistics line from the old file
         szError = '@SEGMENT_MOL: Error reading next atom data'
         read(iUOld,'(a)',err=666,iostat=iErr) szOldBuf(1)
c     end do while (szOldBuf(1)(1:6).ne.'FINISH')
      end do

c   o mark the end of the file
      szError = '@SEGMENT_MOL: Error finalizing new MOL file'
      write(iUNew,'(a)',err=666,iostat=iErr) 'FINISH'

      close(iUOld)
      close(iUNew)
      return
 666  print '(a)', szError
      return
c     end subroutine segment_mol
      end

