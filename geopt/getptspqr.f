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
      SUBROUTINE GETPTSPQR(SCRATCH, EHESS, EIGVH, DMR1, DMR2, HVALUE, 
     &                     EPS, NX, NDIM)
C
C Do a simple binary root seacrh for F(u) = hvalue in the range 
C of U_UPPER < u < U_LOWER. The EPS is the tolerance. 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C

      PARAMETER (MAXTRY = 40, LUOUT = 6)
      LOGICAL UPDATE, NWTONRPSON, NEGEVAL, GENSTPSZ
C
      DIMENSION SCRATCH(NX*NX), EHESS(NDIM, NDIM), DMR1(NDIM, NDIM), 
     &          DMR2(NDIM, NDIM), EIGVH(NDIM, NDIM)
C
      DATA ZERO /0.0D0/, PTO5 /0.5d0/
C
      U_UPPER = 1.0d0
      U_LOWER = 0.0d0
      UNEW    = 0.0d0
      CURVTRE = 0.0d0
      UPDATE  = .FALSE.
      NWTONRPSON = .FALSE.
      GENSTPSZ   = .FALSE.
C
      IF (EHESS(1, 1) .LT. ZERO) THEN
C
         NEGEVAL = .TRUE.
         CALL QSDLINE(SCRATCH, EHESS, EIGVH, DMR1, DMR2, UNEW,
     &                U_UPPER, U_LOWER, HVALUE, CURVTRE, NDIM, NX,
     &                UPDATE, NWTONRPSON, NEGEVAL, GENSTPSZ)

         NEGEVAL = .FALSE. 
C
      ENDIF
C
 10   CONTINUE
C
      UNEW   = (U_UPPER + U_LOWER)/2.0d0
C
      CALL QSDLINE(SCRATCH, EHESS, EIGVH, DMR1, DMR2, UNEW, U_UPPER, 
     &             U_LOWER, HVALUE, CURVTRE, NDIM, NX, UPDATE, 
     &             NWTONRPSON, NEGEVAL, GENSTPSZ)
C
      IF (NWTONRPSON) RETURN
C
      CALL VADD(SCRATCH(1 + 3*NDIM), SCRATCH(1 + 5*NDIM), 
     &          SCRATCH(1 + 4*NDIM), NDIM, -1.00d0)
C
      dtmp = xdot(NDIM,SCRATCH(1 + 3*NDIM),1,SCRATCH(1 + 3*NDIM),1)
      DELTAX = DSQRT(dtmp)
C     
      IF (DABS(DELTAX) .GT. EPS) GO TO 10

C Now calculate the new point!
C
      UPDATE = .TRUE.
      UOLD   = (U_UPPER + U_LOWER)/2.0d0
C
      CALL QSDLINE(SCRATCH, EHESS, EIGVH, DMR1, DMR2, UOLD, U_UPPER,
     &             U_LOWER, HVALUE, CURVTRE, NDIM, NX, UPDATE,
     &             NWTONRPSON, NEGEVAL, GENSTPSZ)
C
      RETURN
      END
