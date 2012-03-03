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
c double  CARTHESS(3*NATOM) : Cartesian coordinate Hessian
c double  DIPDER(9*NATOM)   : Symmetry  coordinate dipole derivative
c double  DIPDFUL(9*NATOM)  : Cartesian coordinate dipole derivative
c double  POLDER(27*NATOM)  : Symmetry  coordinate polarizability derivative
c double  POLDFUL(27*NATOM) : Cartesian coordinate polarizability derivative
c double  DSCR(NDSCR)       : (scr) double scratch

c RECORDS
c get 'NUMPOINT'
c get 'ENGPOINT'
c get 'REFENERG'
c get 'GRDPOINT'
c get 'DIPPOINT'
c get 'POLPOINT'
c get 'OPERSREF'
c get 'NUMVIBRT'
c get 'INVPSMAT'
c get 'NPTIRREP'
c get TYPE//'SYQT'
c get TYPE//'DEGN'
c get TYPE//'LABL'
c put 'HESSIANM'

      SUBROUTINE SETFCM(NATOM,NIRREP,IORDGP,TYPE,
     &                  LABEL,ISYMIRR,INVOP,SYOP,
     &                  SYMHESS,CARTHESS,
     &                  DIPDER,DIPDFUL,
     &                  POLDER,POLDFUL,
     &                  DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 TYPE
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ISYMIRR(3*NATOM),INVOP(3*NATOM),SYOP(9*IORDGP)
      DIMENSION SYMHESS(9*NATOM*NATOM),CARTHESS(9*NATOM*NATOM),
     &          DIPDER(3*3*NATOM),DIPDFUL(3*3*NATOM),
     &          POLDER(3*NATOM*9),POLDFUL(9*3*NATOM)
      double precision dscr(ndscr)

      DIMENSION idegen(100),nptsirr(100)
      LOGICAL PRINTQ, bTmp

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
#include "control.com"

      PRINTQ=(IFLAGS(1).GT.10)
      STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
      NSIZE=3*NATOM
      lFREE=1

      CALL IGETREC(20,'JOBARC','NUMPOINT',1,NPOINT)
      IF (ENERONLY) THEN
         lENGPNT = lFREE
         lFREE   = lFREE + NPOINT
      ELSE
         lGRDPNT = lFREE
         lFREE   = lFREE + NSIZE*NPOINT
         lDIPPNT = lFREE
         lFREE   = lFREE + 3*NPOINT
         lPOLPNT = lFREE
         lFREE   = lFREE + 9*NPOINT
      END IF
      NDSCRLFT = NDSCR+1-lFREE
      IF (NDSCRLFT.LT.0) THEN
         print *, '@SETFCM: Insufficient memory.'
         print *, '         need ',-ndscrlft*ifltln,' more bytes'
         call aces_exit(1)
      END IF
      IF (ENERONLY) THEN
         CALL DGETREC(20,'JOBARC','ENGPOINT',NPOINT,DSCR(lENGPNT))
         Print*, "Energy point"
         Print*, (DSCR(I), I=1, NPOINT)

#ifdef _ASSERT
         do i = 0, npoint-1
            if (dscr(lENGPNT+i).eq.0.d0) then
c            o if any energy is exactly 0, then ACES did not do all points
               print *, '@SETFCM: Assertion failed.'
               print *, '         Energy of point ',1+i,' is 0. a.u.'
               call aces_exit(1)
            end if
         end do
#endif
         CALL DGETREC(20,'JOBARC','REFENERG',1,E0)
         Print*, E0
      ELSE
         CALL DGETREC(20,'JOBARC','GRDPOINT',
     &               NSIZE*NPOINT,DSCR(lGRDPNT))
#ifdef _ASSERT
         do i = 0, nsize*npoint-1, nsize
            bTmp = .true.
            do j = 0, nsize-1
               if (dscr(lGRDPNT+i+j).ne.0.d0) bTmp = .false.
            end do
            if (bTmp) then
c            o if any gradient is exactly 0, then something is PROBABLY wrong
               print *, '@SETFCM: Assertion failed.'
               print *, '         Gradient of point ',1+(i/nsize),
     &                  ' is 0. a.u.'
               call aces_exit(1)
            end if
         end do
#endif
         CALL DGETREC(20,'JOBARC','DIPPOINT',
     &               3*NPOINT,DSCR(lDIPPNT))
         CALL DGETREC(20,'JOBARC','POLPOINT',
     &               9*NPOINT,DSCR(lPOLPNT))
         CALL DGETREC(20,'JOBARC','OPERSREF',IORDGP*9,SYOP)
      END IF

      CALL IGETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      CALL IGETREC(20,'JOBARC','INVPSMAT',NMODE,INVOP)
      CALL IGETREC(20,'JOBARC','NPTIRREP',NIRREP,NPTSIRR)
      CALL IGETREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISYMIRR)
      CALL IGETREC(20,'JOBARC',TYPE//'DEGN',NIRREP,IDEGEN)
      IF (PRINTQ) THEN
      CALL DGETREC(20,'JOBARC',TYPE//'LABL',NIRREP,LABEL)
      END IF

C-----------------------------------------------------------------------
CJDW 6/21/96. Initialize other arguments. Try to make sure intensities
C             come out zero in ENERONLY jobs.
CJDW/AP Jul/98 extended to calculate Raman Intensities.
C-----------------------------------------------------------------------
      CALL ZERO(SYMHESS,9*NATOM*NATOM)
      CALL ZERO(DIPDER, 9*NATOM      )
      CALL ZERO(POLDER, 9*NATOM*3    )
C-----------------------------------------------------------------------

      IFIRST=0
      IFIRSTP=0
      IFIRSTD=0
      IFIRSTG=0
      IOFF=1
      IOFFD=1
      IOFFP=1
      IOFFINV=1
      NLEFT=NMODE

c   o loop over symmetry blocks
      DO IRREP=1,NIRREP

c   o find first occurance of this irrep
      ILOC=ISRCHEQ(NMODE,ISYMIRR,1,IRREP)
      IF (ILOC.NE.NMODE+1) THEN
         ILAST=ISRCHNE(NLEFT,ISYMIRR(ILOC),1,IRREP)
         NVIBSYM=ILAST-1
         NVIBUNQ=NVIBSYM/IDEGEN(IRREP)
         IF (PRINTQ) THEN
            WRITE(6,2000)LABEL(IRREP),IDEGEN(IRREP),NVIBSYM
2000        FORMAT(T3,' Symmetry : ',A,' Degeneracy : ',I1,
     &             ' Unique symmetry coordinates : ',I3)
         END IF

C USE FIRST NVIBUNQ VECTORS FOR DEGENERATE REPS SINCE THEY ARE SORTED
C ACCORDING TO SUBGROUP IRREPS.  OTHERS ARE REDUNDANT.

         IF (ENERONLY) THEN
#ifdef _DEBUG_LVL0
            Print*, "Entering ener2fcm"
#endif
            CALL ENER2FCM(NVIBUNQ,DSCR(lENGPNT+IFIRST),SYMHESS(IOFF),
     &                    INVOP(IOFFINV),STPSIZ,E0,DSCR(lFREE),NDSCRLFT)
            IFIRST=IFIRST+NPTSIRR(IRREP)
         ELSE
C 08/16, Extensions to do Raman Intensities, Ajith and John.
#ifdef _DEBUG_LVL0
            Print*, "Entering the grad2fcm"
#endif
            CALL GRAD2FCM(NATOM,NVIBUNQ,
     &                    DSCR(lGRDPNT+IFIRSTG),SYMHESS(IOFF),
     &                    DSCR(lDIPPNT+IFIRSTD),DIPDER(IOFFD),
     &                    DSCR(lPOLPNT+IFIRSTP),POLDER(IOFFP),
     &                    INVOP(IOFFINV),SYOP,STPSIZ,
     &                    DSCR(lFREE),NDSCRLFT)
            IOFFD=IOFFD+3*NVIBUNQ
            IOFFP=IOFFP+9*NVIBUNQ
            IFIRSTG=IFIRSTG+NPTSIRR(IRREP)*NSIZE+NVIBSYM
            IFIRSTD=IFIRSTD+NPTSIRR(IRREP)*3
            IFIRSTP=IFIRSTP+NPTSIRR(IRREP)*9
         END IF
         IOFF=IOFF+NVIBUNQ*NVIBUNQ
         IOFFINV=IOFFINV+NVIBSYM
         NLEFT=NLEFT-NVIBUNQ

      END IF
      END DO

c   o reset the scratch (l*PNT is not used anymore)
      lFREE = 1
      NDSCRLFT = NDSCR

c   o transform and write the Hessian to JOBARC
      CALL TRNFCM(NMODE,NATOM,NIRREP,
     &            SYMHESS,CARTHESS,
     &            ISYMIRR,IDEGEN,TYPE,DSCR(lFREE),NDSCRLFT)
      CALL DPUTREC(20,'JOBARC','HESSIANM',NSIZE*NSIZE,CARTHESS)

C Modifications for Raman Intensities Ajith and John 08/98
      CALL TRNDIP(NMODE,NATOM,NIRREP,
     &            DIPDER,DIPDFUL,POLDER,POLDFUL,
     &            ISYMIRR,IDEGEN,TYPE,DSCR(lFREE),NDSCRLFT)

      RETURN
      END

