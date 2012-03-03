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

      SUBROUTINE MSP_PSB(V, H, SCRATCH, STEP, TBT, TOLRNCE,
     &                   POWEL_DONE, MSPORPSB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL POWEL_DONE
      CHARACTER*6 MSPORPSB

#include "io_units.par"
#include "fnamelen.par"
#include "jodaflags.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION V(NXM6),H(NXM6*NXM6),SCRATCH(NX*NX),TBT(3*NXM6*NXM6),
     $          STEP(NXM6)
      DATA ONE /1.00D0/
C
C Murtagh-Sargent update. Note that there are conflicting formulas
C MS update. The conflicts are about the scaling factors and I
C doubt that they can cause major problems.
C
      N2 = NXM6 + 1
      N3 = 2*NXM6+1
      N4 = 3*NXM6+1

C SCRATCH(N2) = HDX where H is the Hessian
C SCRATCH(N3) = DG (gradient difference)
C SCRATCH(N4) = (DG - HDX)

      DDXDX  = DDOT(NXM6, STEP, 1, STEP, 1)
      DDXDX2 = DDXDX**2
      DHDXDG = DDOT(NXM6, SCRATCH(N4), 1, SCRATCH(N4), 1)
      DHDXDX = DDOT(NXM6, SCRATCH(N4), 1, STEP, 1)

      IF (MSPORPSB .EQ. "BOFILL") THEN
         IF (DABS(DHDXDG).GT.TOLRNCE.AND.
     &       DABS(DHDXDX).GT.TOLRNCE) THEN
            COEF    = (DHDXDX**2)/(DDXDX*DHDXDG)
            PHI_INT = ONE/DHDXDX
         ELSE
            COEF    = ONE
            PHI_INT = ONE
         END IF
      ELSE IF (MSPORPSB .EQ. "MSP   ") THEN
         PHI_INT =  ONE/DHDXDX
      ELSE IF (MSPORPSB .EQ. "PSB   ") THEN
         PHI_INT = DHDXDX/DDXDX2
      END IF

      N2 = NXM6*NXM6+1
      N3 = NXM6*NXM6+N2
C
C Form (DG - HDX)o(DG - HDX) and put it in SCRATCH(N2)
C
      CALL MODMATMUL(SCRATCH(N2), SCRATCH(N4), SCRATCH(N4), NXM6, 1,
     &               NXM6, NXM6, 1, NXM6)

      IF (MSPORPSB .EQ. "BOFILL") THEN
         CALL ZERO(TBT(N3), NXM6*NXM6)
         CALL XDAXPY(NXM6*NXM6, PHI_INT, SCRATCH(N2), 1, TBT(N3), 1)
         CALL XDAXPY(NXM6*NXM6, COEF, TBT(N3), 1, H, 1)
      ELSE IF (MSPORPSB .EQ. "MSP   ") THEN
         CALL XDAXPY(NXM6*NXM6, PHI_INT, SCRATCH(N2), 1, H,  1)
      END IF

      IF (.NOT. POWEL_DONE .AND. (MSPORPSB .EQ. "PSB   ")) THEN
C        Form DXoDX in TBT(1)
         CALL MODMATMUL(TBT(1),STEP,STEP,NXM6,1,NXM6,NXM6,1,
     &                  NXM6)
C        Form (DG - HDX)DX in TBT(N2)
         CALL MODMATMUL(TBT(N2),SCRATCH(N4),STEP,NXM6,1,NXM6,
     &                  NXM6,1,NXM6)
C        TBT(N3) = TBT(N2)^{trnsp}
         CALL MTRANSP(TBT(N2),TBT(N3),NXM6,NXM6,NXM6,NXM6)
C        TBT(N3) = (DG - HDX)DX + ((DG - HDX)DX)^{trnsp}
         CALL VADD(TBT(N3),TBT(N3),TBT(N2),NXM6*NXM6,1.D0)
         CALL XDSCAL(NXM6*NXM6, ONE/DDXDX, TBT(N3), 1)
         CALL XDAXPY(NXM6*NXM6, -PHI_INT, TBT(1), 1, TBT(N3), 1)
         CALL XDAXPY(NXM6*NXM6, ONE, TBT(N3), 1, H, 1)
      END IF

      IF (IFLAGS(1).GE.10) THEN
         IF (MSPORPSB .EQ. "BOFILL") THEN
             Print*, "-----The Bofill updated Hessian-----"
         ELSE IF (MSPORPSB .EQ. "PSB   ")  THEN
             Print*, "-----The PSB updated Hessian-----"
         ELSE IF (MSPORPSB .EQ. "MSP   ") THEN
             Print*, "-----The MSP updated Hessian-----"
         ENDIF
         CALL HESSOUT(H,NXM6,NXM6,0)
      END IF

      RETURN
      END

