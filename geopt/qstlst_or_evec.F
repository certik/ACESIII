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

      SUBROUTINE QSTLST_OR_EVEC(QSTLST_TANGENT, GRDMOD, HESMOD, DIAGHES,
     &                          HES, SCRATCH, IMODE, QSTLST_CLIMB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL QSTLST_CLIMB
      PARAMETER (THRESHOLD = 5.0D-2, WGHT= 1.0D-3)

#include "mxatms.par"
#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DATA HALF /0.50D0/

      DIMENSION HESMOD(NOPT, NOPT), GRDMOD(NOPT), QSTLST_TANGENT(NOPT),
     &          SCRATCH(NX*NX), DIAGHES(NOPT, NOPT),
     &          HES(NOPT, NOPT)

      CALL ZERO(SCRATCH, NOPT)
c      Print*, (QSTLST_TANGENT(I), I=1,NOPT)
      CALL XGEMM('N','N', NOPT, 1, NOPT, 1.0D0, HES, NOPT,
     &            QSTLST_TANGENT, NOPT, 0.0D0, SCRATCH, NOPT)

      EIGVALUE = DDOT(NOPT, QSTLST_TANGENT, 1, SCRATCH, 1)
      GRDVALUE = DDOT(NOPT, QSTLST_TANGENT, 1, GRDMOD,  1)

c      Write(6,*) "Norm of the projected grad and Hes", EIGVALUE,
c     &            GRDVALUE
      EST_STEP = (EIGVALUE + DSQRT(EIGVALUE**2 + 4.0D0*(GRDVALUE)
     &            **2))*HALF
c      Write(6,*) "Debug-Info1: The estimated step size", EST_STEP
      EST_STEP =  -GRDVALUE/(EIGVALUE - EST_STEP)

c      Write(6,*) "Debug-Info2: The estimated step size", EST_STEP
C
C Let's check that the estimated step size is greater than 0.05 au.
C If so, then follow the QST or LST tangent (I am interpreting Schlegel's
C "estimated displacement along the tangent vector is greater than the 0.05".
C At this point I am not sure whether he meant the largest absolute
C displacement. Also cap the QST/LST to first four cycle (also rec. by Schlegel)
C
      QSTLST_CLIMB = ((DABS(EST_STEP).GE.THRESHOLD.AND.NCYCLE.LE.4).OR.
     &                NCYCLE.LE.2)
c      Write(6,*) "The climbing Phase if true", QSTLST_CLIMB
C
C Compute the gradient along the tangent vector, copy the tangent
C vector to the eigenvector matrix (to the first eigenvector) and the
C eigenvalue (EIGVALUE) to the first eigenvalue and set the IMODE
C to 1.
C
CSSS      IF (QSTLST_CLIMB) THEN
CSSS         WEIGHT = ONE - WGHT
CSSS         CALL MODFY_HESSIAN(DIAGHES, HESMOD, HES, QSTLST_TANGENT,
CSSS     &                      SCRATCH, EIGVALUE, WEIGHT, NOPT)
CSSS      END IF
C
C This is just to be safe in subsequent steps.
C
      CALL ZERO(SCRATCH,NOPT*NOPT)

      RETURN
      END

