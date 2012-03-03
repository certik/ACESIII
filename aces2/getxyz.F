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

C Read Cartesian coordinates from Z-matrix

      SUBROUTINE GETXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#include "mxatms.par"
#include "io_units.par"
#include "linelen.par"
#include "flags.h"
      DOUBLE PRECISION BTOA
      PARAMETER (BTOA=0.529177249d0)

C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     1     PARNAM(3*MXATMS)
#include "coord_dup.com"
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
#include "machsp.com"

C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol

      double precision xyz(3,mxatms)
      integer izl(2,7), iStruct
      CHARACTER*(linelen) ZLINE
      logical btmp, bStruct, bAtom

c ----------------------------------------------------------------------

c   o scaling factor for Angstoms to Bohr
      FACTOR = 1.d0 / BTOA

c   o open up the main input file
      open(LuZ,FILE=ZFil,FORM='FORMATTED',STATUS='OLD')
      rewind LuZ

c   o skip the header (zline=TITLE on end do)
      btmp = .true.
      do while (btmp)
         read(luz,'(a)') zline
         call parsez(zline,izl)
         i = izl(1,1)
         btmp = (i.eq.0).or.(zline(i:i).eq.'%')
      end do
      call putcrec(1,' ','TITLE',80,zline(1:80))

c   o load the first structure (reactant) into Q and ZSYM
      ndx = 0
      bAtom = .true.
      do while (bAtom)
         read(luz,'(a)') zline
         call parsez(zline,izl)
         if (izl(1,1).ne.0) then
c         o increment atom counter
            ndx = ndx + 1
            if (ndx.gt.mxatms) then
               print *, '@GETXYZ: Primary structure exceeds the ',
     &                  'maximum number of atoms (',mxatms,')'
               call errex
            end if
c         o load the symbol
            if (izl(2,1).gt.izl(1,1)+1) then
               print *, '@GETXYZ: Symbol token is invalid'
               print *, '         atom ',ndx,' = "',
     &                  zline(izl(1,1):izl(2,1)),'"'
               call errex
            end if
            zsym(ndx)=zline(izl(1,1):izl(2,1))
c         o load the coordinates
            if (izl(1,2).eq.0.or.izl(1,3).eq.0.or.
     &          izl(1,4).eq.0.or.izl(1,5).ne.0    ) then
               print *, '@GETXYZ: Error in XYZ coordinates'
               print *, '         atom ',ndx,': ',zline
               call errex
            end if
            read(zline(izl(1,2):izl(2,4)),*) (q(i,ndx),i=1,3)
         else
            bAtom = .false.
         end if
c     end do while (bAtom)
      end do

c   o process the reactant as the main structure
      if (ndx.eq.0) then
         print *, '@GETXYZ: Missing XYZ coordinates'
         call errex
      end if
      natoms = ndx
      nx = 3*natoms
      if (iflags(h_IOPPAR_units).eq.0) call dscal(nx,factor,q,1)
      call pertable

c   o flush remaining atom slots
      if (natoms.lt.mxatms) then
         i = 3*(mxatms-natoms)
         call zero(q(1,natoms+1),i)
         do i = natoms+1, mxatms
            zsym(i) = ' '
         end do
      end if

c   o read additional structures into xyz until we hit a namelist
c     (iStruct corresponds to the structure we expect to read)
      iStruct = 1
      bStruct = .true.
      do while (bStruct.and.iStruct.lt.3)
         iStruct = iStruct + 1
         ndx = 0
         bAtom = .true.
         do while (bAtom)
c         o only read the first line of the secondary structure
            if (ndx.ne.0.or.iStruct.eq.2) then
               read(luz,'(a)') zline
               call parsez(zline,izl)
            end if
            i = izl(1,1)
            if (i.ne.0.and.(ndx.ne.0.or.zline(i:i).ne.'*')) then
c            o increment atom counter
               ndx = ndx + 1
               if (ndx.gt.natoms) then
                  print *, '@GETXYZ: Structure ',iStruct,' has more ',
     &                     'atoms than the primary structure.'
                  print *, '         atom ',ndx,': ',zline
                  call errex
               end if
c            o make sure the symbol is the same
               if (zsym(ndx).ne.zline(izl(1,1):izl(2,1))) then
                  print *, '@GETXYZ: Atoms do not match'
                  print *, '         ',zsym(ndx),' != "',
     &                     zline(izl(1,1):izl(2,1)),'"'
                  call errex
               end if
c            o load the coordinates
               if (izl(1,2).eq.0.or.izl(1,3).eq.0.or.
     &             izl(1,4).eq.0.or.izl(1,5).ne.0    ) then
                  print *, '@GETXYZ: Error in XYZ coordinates'
                  print *, '         structure ',iStruct
                  print *, '         atom ',ndx,': ',zline
                  print *,'izl ',(izl(1,i),i=1,5)
                  call errex
               end if
               read(zline(izl(1,2):izl(2,4)),*) (xyz(i,ndx),i=1,3)
            else
               bAtom = .false.
               if (ndx.eq.0) then
c               o we hit a namelist or another blank line
                  bStruct = .false.
               else
c               o compare this structure to the primary structure
                  if (ndx.ne.natoms) then
                     print *, '@GETXYZ: Structure ',iStruct,' has ',
     &                        'fewer atoms than structure 1'
                     call errex
                  end if
               end if
            end if
c        end do while (bAtom)
         end do
         if (bStruct) then
            if (iStruct.eq.2) then
c            o read the next line and see if this structure was the TS or PR
               read(luz,'(a)') zline
               call parsez(zline,izl)
               i = izl(1,1)
               btmp = (i.eq.0.or.zline(i:i).eq.'*')
            else
               btmp = .true.
            end if
            if (iprnt.gt.1
#ifdef _DEBUG
     &          .or..true.
#endif
     &         ) then
               if (btmp) then
                  print *, 'Product structure:'
               else
                  print *, 'Transition-state structure:'
               end if
               do ndx = 1, natoms
                  print *, zsym(ndx),(xyz(i,ndx),i=1,3)
               end do
c           end if (print structures)
            end if
            if (iflags(h_IOPPAR_units).eq.0) 
     *              call dscal(nx,factor,xyz,1)
            if (btmp) then
               call dputrec(1,'JOBARC','PRSTRUCT',nx,xyz)
            else
c            o move Q to RXSTRUCT since we are optimizing the TS structure
               call dputrec(1,'JOBARC','RXSTRUCT',nx,q)
               call dcopy(nx,xyz,1,q,1)
            end if
c        end if (bStruct)
         end if
c     end do while (bStruct)
      end do

      call geomout(natoms, q, zsym, iatnum)
      RETURN
      END

