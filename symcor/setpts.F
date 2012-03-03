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
C     we need to make appropriate allowances.
C-----------------------------------------------------------------------

c INPUT
c integer NATOM
c integer NIRREP
c char*4  TYPE
c integer NDSCR

c OUTPUT
c integer NENER
c char*8  LABEL(NIRREP)
c integer ISYMIRR(3*NATOM)
c integer IPTTYPE(9*NATOM*NATOM)
c integer INVOP(3*NATOM)
c double  POINTS(27*NATOM*NATOM*NATOM)
c double  DSCR(NDSCR)

c RECORDS
c get 'NFDIRREP'
c get 'FDIRREP '
c get 'COORD   '
c get 'ATOMMASS'
c get TYPE//'SYMQ'
c get TYPE//'SYQT'
c get TYPE//'DEGN'
c get TYPE//'LABL'
c get 'INVPSMAT'
c get 'NUMVIBRT'
c put 'NUMPOINT'
c put 'FDCALCTP'
c put 'FDCOORDS'
c put 'NPTIRREP'
c put 'ENGPOINT'
c put 'GRDPOINT'
c put 'DIPPOINT'
c put 'POLPOINT'

      SUBROUTINE SETPTS(NATOM,NIRREP,TYPE,
     &                  NENER,LABEL,ISYMIRR,IPTTYPE,INVOP,
     &                  POINTS,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 TYPE
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ISYMIRR(3*NATOM),IPTTYPE(9*NATOM*NATOM),INVOP(3*NATOM)
      DIMENSION POINTS(27*NATOM*NATOM*NATOM)
      DIMENSION DSCR(NDSCR)

      INTEGER SKIPIR
      DIMENSION IJUNK(20),SKIPIR(20),idegen(100),nptirr(100)
      LOGICAL PRINTQ

      COMMON /FLAGS/  IFLAGS(100)
#include "control.com"
#include "machsp.com"

      DATA TOL /1.D-8/

      PRINTQ=(IFLAGS(1).GE.10)
      NSIZE=3*NATOM

      irank = 0   ! simulate single-processor execution under ACESII.
      nprocs = 1

      if (gmtryopt.or.gradonly.or.nprocs.eq.1) then
         iRoot = 0
      else
c      o the root process is doing the reference geometry in a vib freq calc
c        with numerical gradients
         iRoot = 1
      end if

c   o check for restrictions on FD irreps
      CALL IGETREC(-1,'JOBARC','NFDIRREP',1,ICOUNT)
      IF (ICOUNT.NE.0) THEN
         CALL IGETREC(-1,'JOBARC','FDIRREP ',ICOUNT,IJUNK)
         WRITE(6,501)
501      FORMAT(T3,'@SETPTS: FCM evaluation limited ',
     &             'to the following symmetries :')
         WRITE(6,'((12I5))')(IJUNK(I),I=1,ICOUNT)
         DO I=1,20
            SKIPIR(I)=1
         END DO
         DO I=1,ICOUNT
            SKIPIR(IJUNK(I))=0
         END DO
      ELSE
         CALL IZERO(SKIPIR,20)
      END IF

      STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
      IF (PRINTQ) THEN
         WRITE(6,500)STPSIZ
500      FORMAT(T3,'Step size will be ',F8.5,' amu**(1/2) * bohr.')
      END IF
#ifdef _DEBUG_LVL0
      Print*, "The step size in symcor:", STPSIZ
#endif 
      lFREE = 1

      lSYMQ = lFREE
      lFREE = lFREE + NSIZE*NSIZE
      lCOOR = lFREE
      lFREE = lFREE + NSIZE
      lMASS = lFREE
      lFREE = lFREE + NATOM
      NDSCRLFT = NDSCR+1-lFREE
      IF (NDSCRLFT.LT.0) THEN
         print *, '@SETPTS: Insufficient memory.'
         print *, '         need ',-ndscrlft*ifltln,' more bytes'
         call aces_exit(1)
      END IF
      CALL DGETREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE,
     &                                                DSCR(lSYMQ))
      CALL DGETREC(20,'JOBARC','COORD   ',NSIZE,DSCR(lCOOR))
#ifdef _DEBUG_LVL0
      Print*, "The Cartesian coordinates read in septs:"
      Print*, (Dscr(I), I=lCOOR, NSIZE)
#endif      
      CALL DGETREC(20,'JOBARC','ATOMMASS',NATOM,DSCR(lMASS))

      CALL IGETREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISYMIRR)
      CALL IGETREC(20,'JOBARC',TYPE//'DEGN',NIRREP,IDEGEN)
      CALL IGETREC(20,'JOBARC','INVPSMAT',3*NATOM,INVOP)
      IF (PRINTQ) THEN
      CALL DGETREC(20,'JOBARC',TYPE//'LABL',NIRREP,LABEL)
      END IF

      CALL IGETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      NLEFT=NMODE

      CALL ZERO(POINTS,27*NATOM*NATOM*NATOM)
      CALL IZERO(NPTIRR,NIRREP)

c   o set up vector of reciprocal square roots of atomic masses
      DO IOFF=0,NATOM-1
         X=SQRT(DSCR(lMASS+IOFF))
         IF (X.LT.TOL) THEN
            DSCR(lMASS+IOFF)=0.d0
         ELSE
            DSCR(lMASS+IOFF)=1.d0/X
         END IF
      END DO

c   o loop over irreducible representations and process first occurances
      NENER=0
      NGRAD=0
      NPTOTX=0
      NENERX=0
      NGRADX=0
      IOFF=0
      IPOS2=1
      IPOS3=1
      DO IRREP=1,NIRREP
         IFIRST=ISRCHEQ(NMODE,ISYMIRR,1,IRREP)
      IF (IFIRST.NE.NMODE+1) THEN

         ILAST=ISRCHNE(NLEFT,ISYMIRR(IFIRST),1,IRREP)
         NVIBSYM=ILAST-1
         NVIBUNQ=NVIBSYM/IDEGEN(IRREP)

         IF (PRINTQ) THEN
            WRITE(6,2000)LABEL(IRREP),IDEGEN(IRREP),NVIBUNQ
2000        FORMAT(T3,' Symmetry : ',A,' Degeneracy : ',I1,
     &             ' Unique symmetry coordinates : ',I3)
         END IF

         IPOS=NSIZE*(IFIRST-1)
#ifdef _DEBUG_LVL0
         Print*, "In setpts, energy only?", ENERONLY 
#endif
         IF (ENERONLY) THEN
            CALL DOENER(NATOM,NVIBUNQ,
     &                  DSCR(lSYMQ+IPOS),DSCR(lCOOR),STPSIZ,DSCR(lMASS),
     &                  POINTS(IPOS2),NPOINT,
     &                  INVOP(IPOS3),PRINTQ,DSCR(lFREE),NDSCRLFT)
            NENER=NENER+NPOINT
            IF (SKIPIR(IRREP).EQ.0) THEN
               do iR = 0, nProcs-1
                  call paces_batch_stat(iR,nProcs,iRoot,nPoint,iO,nE)
                  do iX = 1+iO, nE+iO
                     iPtType(iOff+iX) = 1+iR
                  end do
                  if (iR.eq.iRank) nEnerX=nEnerX+nE
               end do
               iRoot = mod(iRoot+nPoint,nProcs)
               nPTotX=nPTotX+nPoint
            ELSE
               DO IX=1,NPOINT
                  IPTTYPE(IOFF+IX)=0
               END DO
            END IF
         ELSE
#ifdef _DEBUG_LVL0
         Print*, "In setpts, grad only?, T"
#endif  

            CALL DOGRAD(NATOM,NVIBUNQ,
     &                  DSCR(lSYMQ+IPOS),DSCR(lCOOR),STPSIZ,DSCR(lMASS),
     &                  POINTS(IPOS2),NPOINT,
     &                  INVOP(IPOS3),PRINTQ,DSCR(lFREE),NDSCRLFT)
            NGRAD=NGRAD+NPOINT
            IF (SKIPIR(IRREP).EQ.0) THEN
               do iR = 0, nProcs-1
                  call paces_batch_stat(iR,nProcs,iRoot,nPoint,iO,nE)
                  do iX = 1+iO, nE+iO
                     iPtType(iOff+iX) = 1+iR
                  end do
                  if (iR.eq.iRank) nGradX=nGradX+nE
               end do
               iRoot = mod(iRoot+nPoint,nProcs)
               nPTotX=nPTotX+nPoint
            ELSE
               DO IX=1,NPOINT
                  IPTTYPE(IOFF+IX)=0
               END DO
            END IF
         END IF
         NPTIRR(IRREP)=NPTIRR(IRREP)+NPOINT
         IPOS2=IPOS2+NPOINT*NSIZE
         IPOS3=IPOS3+NVIBSYM
         IOFF=IOFF+NPOINT
         NLEFT=NLEFT-NVIBSYM

      END IF
      END DO

      NTOT=NGRAD+NENER
      NTOTX=NGRADX+NENERX
      WRITE(6,1000)NTOTX,NENERX,NGRADX
1000  FORMAT(T3,' Total number of calculations required      : ',I5,/,
     &       T3,' Number of single-point energy calculations : ',I5,/,
     &       T3,' Number of energy gradient     calculations : ',I5)
      IF (NTOT.GT.NSIZE*NSIZE) THEN
         WRITE(6,*) '@SETPTS: Too many points!'
         CALL ERREX
      END IF
      IF (NPTOTX.EQ.0) THEN
         WRITE(6,*) '@SETPTS: There are no vibrational modes!'
         CALL ERREX
      END IF
      IF (NTOTX.EQ.0) THEN
         if (nprocs.eq.1.or.irank.ne.0.or.gmtryopt.or.gradonly) then
            WRITE(6,*) '@SETPTS: There are no points to calculate!'
            CALL ERREX
         else
            CALL IPUTREC(20,'JOBARC','LASTGEOM',1,1)
         end if
      END IF

c   o write geometries to be used and calculation types to jobarc
      CALL IPUTREC(20,'JOBARC','NUMPOINT',1,NTOT)
      CALL IPUTREC(20,'JOBARC','FDCALCTP',NTOT,IPTTYPE)
      CALL DPUTREC(20,'JOBARC','FDCOORDS',NSIZE*NTOT,POINTS)
      CALL IPUTREC(20,'JOBARC','NPTIRREP',NIRREP,NPTIRR)

c   o write out zero vectors for calculation archiving
      CALL ZERO(POINTS,NTOT*max(NSIZE,9))
      CALL DPUTREC(20,'JOBARC','ENGPOINT', NTOT,POINTS)
      CALL DPUTREC(20,'JOBARC','GRDPOINT',NSIZE*NTOT,POINTS)
      CALL DPUTREC(20,'JOBARC','DIPPOINT', 3*NTOT,POINTS)
      CALL DPUTREC(20,'JOBARC','POLPOINT', 9*NTOT,POINTS)

      RETURN
      END

