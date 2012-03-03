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
      SUBROUTINE HUPDATE(V,H,SCRATCH,STEP,TBT,ITYPE)
C
C     ARGUMENTS PASSED ARE NEW GRADIENT, HESSIAN (TO BE READ AND
C     RETURNED UPDATED), SCRATCH ARRAY- OLD GRADIENT READ INTO FIRST
C     NXM6 POSITIONS, UPDATE TYPE
C     TBT IS JUST A SCRATCH ARRAY FOR MATRIX MULTIPLY OPERATIONS
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "io_units.par"
#include "fnamelen.par"
#include "jodaflags.com"
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION V(NXM6),H(NXM6*NXM6),SCRATCH(NX*NX),TBT(3*NXM6*NXM6),
     $   STEP(NXM6)
      LOGICAL POWEL_DONE
      CHARACTER*6 MSPORPSB
      CHARACTER*10 TRASH
      CHARACTER*(fnamelen) FNAME
      INTEGER    LuArc
      PARAMETER (LuArc  = 77)
      CHARACTER*(*) ArcFil
      PARAMETER    (ArcFil = 'OPTARC ')
      PARAMETER (TOLRNCE = 1.0D-12)
      DATA ONE /1.00D0/
C
      CALL GFNAME(ARCFIL,FNAME,ILENGTH)
      OPEN(UNIT=LUARC,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
      DO 10 I=1,NCYCLE
 10   READ(LUARC)TRASH
C
C     Most of the stuff read here is not used - hence the stupid vars.
C The previous cycle's gradients and Hessians are read into scratch and H.
C The incoming variables V, H, and STEP contain the current gradient, Hessian,
C and step size.
C
      READ(LUARC)TRASH,ICRAP,icycle,NJ,NJN,NJNK,NOPT,NJ1,NJ2,NJ3,NJ4,
     $   NJ5,NJ6,NJ7,NJ8,NJ9,NJA,NJB,NJC,NJD,NJE,NJF,NJG,E,
     &   (SCRATCH(I),I=1,NXM6),
     &   (H(I),I=1,NXM6*NXM6),(STEP(J),J=1,NXM6)
      CLOSE(UNIT=LUARC)
      IF(ITYPE.EQ.0)THEN
         WRITE(LuOut,70) icycle
 70      FORMAT(T3,' Hessian from cycle ',i2,' read and not updated.')
      ELSEIF(ITYPE.EQ.1)THEN
C
C        Powell update
C
         WRITE(LuOut,71) icycle
 71      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Powell update using ',
     &      'last two gradients and previous step.')
C
         PHI = ONE
         CALL POWEL(V, H, SCRATCH, STEP, TBT, PHI)
C
      ELSEIF(ITYPE.EQ.2)THEN
C
C        Broyden-Fletcher-Goldfarb-Shanno update
C
         WRITE(LuOut,72) icycle
 72      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' BFGS update using ',
     &      'last two gradients and previous step.')
C
#ifdef _DEBUG_LVLM1
         Print*, "The BFGS update check"
         Write(6,*) (V(I), I =1, NXM6)
         Write(6,*)
         Write(6,*) (SCRATCH(I), I = 1, NXM6)
         Write(6,*)
         Write(6,*) (STEP(I), I= 1, NXM6)
         Write(6,*)
         CALL OUTPUT(H, 1, NXM6, 1, NXM6, NXM6,NXM6, 1 )
         Write(6,*)
#endif
         CALL BFGS(V, H, SCRATCH, STEP, TBT)

C
      ELSEIF(ITYPE.EQ.3)THEN
C
C        Murtagh-Sargent update
C
         WRITE(LuOut,73) icycle
 73      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Murtagh-Sargent update using ',
     &      'last two gradients and previous step.')
C
         POWEL_DONE = .TRUE.
         MSPORPSB   = "MSP   "
         CALL MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                POWEL_DONE, MSPORPSB)
C
      ELSEIF(ITYPE.EQ.4)THEN
C
C        (1-phi)*Powell + (phi)*Murtagh-Sargent update
C
         WRITE(LuOut,74) icycle
 74      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Bofill mixture of Powell and Murtagh-Sargent',
     &      ' update using last two gradients and previous step.')
C
         N2 = NXM6 + 1
         N3 = 2*NXM6+1
         N4 = 3*NXM6+1

         CALL MODMATMUL(SCRATCH(N2),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
         CALL DCOPY(NXM6, V, 1, SCRATCH(N3), 1)
         CALL XDAXPY(NXM6, -ONE, SCRATCH(1), 1, SCRATCH(N3), 1)
         CALL VADD(SCRATCH(N4),SCRATCH(N3),SCRATCH(N2),NXM6,-1.D0)

         DDXDX  = XDDOT(NXM6, STEP, 1, STEP, 1)
         DHDXDG = XDDOT(NXM6, SCRATCH(N4), 1, SCRATCH(N4), 1)
         DHDXDX = XDDOT(NXM6, SCRATCH(N4), 1, STEP, 1)

         IF (DABS(DHDXDX).GT.TOLRNCE.AND.
     &       DABS(DHDXDG).GT.TOLRNCE) THEN
            PHI = (DHDXDX**2)/(DDXDX*DHDXDG)
         ELSE
            PHI = ONE
         END IF

         CALL POWEL(V, H, SCRATCH, STEP, TBT, ONE-PHI)
         POWEL_DONE = .TRUE.
         MSPORPSB   = "BOFILL"
         CALL MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                POWEL_DONE, MSPORPSB)

      ELSE IF (ITYPE .EQ. 5) THEN
C
C        Powell-symmetric-Broyden update
C
         WRITE(LuOut,75) icycle
 75      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' Powell-Symmetric-Broyden update using ',
     &      'last two gradients and previous step.')
C
         POWEL_DONE = .FALSE.
         MSPORPSB   = "PSB   "
         CALL MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                POWEL_DONE, MSPORPSB)

      ELSE IF (ITYPE .EQ. 6) THEN
C
C  Modified BFGS update

         WRITE(LuOut,76) icycle
 76      FORMAT(T3,' Hessian from cycle ',i2,' read.',/,t3,
     $      ' modified BFGS update using ',
     &      'last two gradients and previous step.')
C
#ifdef _DEBUG_LVL0
         Print*, "The BFGS update check"
         Write(6,*) (V(I), I =1, NXM6)
         Write(6,*) (SCRATCH(I), I = 1, NXM6)
         Write(6,*) (STEP(I), I= 1, NXM6)
         CALL OUTPUT(H, 1, NXM6, 1, NXM6, NXM6,NXM6, 1 )
         Write(6,*)
#endif
CSSS         CALL MODF_BFGS(V, H, SCRATCH, STEP, TBT)
      ENDIF
C
      RETURN
      END
