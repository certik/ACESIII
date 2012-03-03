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
      SUBROUTINE CHKECK(NRX,HBUF,ECKART)
C
C SUBROUTINE CHECKS NORMAL MODES AGAINST ECKART CONDITIONS.  PUTS
C   VALUES OF ROTATIONAL AND TRANSLATIONAL PSEUDOMOMENTA INTO A VECTOR
C   AND PASSES THIS BACK TO THE CALLING PROGRAM.  PRINTING IS AN OPTION.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "coord.com"
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
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
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
 
 
      Common /Orient/ Orient(3,3)
#include "io_units.par"
      DIMENSION ECKART(2*NRX),HBUF(NRX,NRX),ZAT(3),ZA(3)
C
C CHECK EIGENVECTORS AGAINST ECKART CONDITIONS - CAN DETERMINE IF
C  IT IS A TRANSLATION, ROTATION OR VIBRATION.  USE TOLERANCES FOR
C  LINEAR AND ANGULAR MOMENTUM.  PUT VALUE INTO ECKART VECTOR
C  (TOP HALF USED FOR ROTATIONS, BOTTOM FOR TRANSLATIONS).
C
C
C FIRST, MAKE EIGENVECTOR DIMENSION = MASS * DISPLACEMENT SO THAT ROTATI
C  ECKART CONDITION CAN BE TESTED LATER.
C
 3456 DO 32 I=1,NRX
      DO 33 J=1,NRX
      HBuf(J,I)=HBuf(J,I)*DSQRT(ATMASS(1+(J-1)/3))
   33 CONTINUE
      CALL NORMAL(HBUF(1,I),NRX)
   32 CONTINUE
C
C TEST ECKART CONDITIONS.
C
      CALL ZERO(ECKART,2*NRX)
      DO 35 I=1,NRX
      CALL ZERO(ZAT,3)
C
C 1. CHECK ANGULAR MOMENTUM CRITERIA
C
       DO 36 J=1,NRX-2,3
        CALL CROSS(HBUF(J,I),Q(J),ZA,0)
        CALL VADD(ZAT,ZAT,ZA,3,1.D0)
   36  CONTINUE
       CALL ZERO(ZA,3)
       ECKART(I)=DIST(ZAT,ZA)
C
C 2. CHECK LINEAR MOMENTUM CRITERIA
C
C
        CALL ZERO(ZAT,3)
        DO 39 J=1,NRX-2,3
        CALL VADD(ZAT,ZAT,HBUF(J,I),3,1.D0)
   39  CONTINUE
       ECKART(NRX+I)=DIST(ZAT,ZA)
   38  CONTINUE
   35  CONTINUE
C
C DUMP ROTATIONAL AND VIBRATIONAL VALUES IF CALLED FOR.
C
      IF(IPRNT.GE.8)THEN
       WRITE(LUOUT,5700)
 5700  FORMAT(T3,' Linear and Angular "momenta" corresponding ',
     &'to Hessian Eigenvectors ',/,T10,'Mode',T30,'P',T50,'L')
       WRITE(LUOUT,5701)(I,ECKART(NRX+I),ECKART(I),I=1,NRX)
 5701  FORMAT(T10,I3,T23,F16.12,T43,F16.12)
      ENDIF
      RETURN
      END
