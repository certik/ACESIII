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
      SUBROUTINE INERTIA(IT, NOSILENT)
C
C BUILD INERTIA TENSOR
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "coord.com"
C
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
 
 
      Parameter (LuOut = 6)
      DOUBLE PRECISION IT(3,3),SCRATCH(9)
      LOGICAL NOSILENT
C
      IND(I,J)=3*(I-1)+J
C
      CALL ZERO(SCRATCH,9)
      CALL ZERO(IT,9)
      DO 100 I=1,3
         DO 100 K=1,NATOMS
            X1=DIST(Q(3*K-2),SCRATCH(1))
            X2=Q(3*(K-1)+I)
 100        IT(I,I)=ATMASS(K)*(X1*X1-X2*X2)+IT(I,I)
c     $         (DIST(Q(3*K-2),SCRATCH(1))**2-Q(3*(K-1)+I)**2)+IT(I,I)
      DO 200 I=1,3
         DO 300 J=1,3
            IF(I.EQ.J)GOTO 300
            DO 400 K=1,NATOMS
               IT(I,J)=-ATMASS(K)*(Q(IND(K,I))*Q(IND(K,J)))+IT(I,J)
 400        Continue
 300     Continue
 200  CONTINUE
      IF(IPRNT.GE.4 .AND. NOSILENT)WRITE(LuOut,*)' Inertia tensor'
      IF(IPRNT .GE. 4 .AND. NOSILENT) WRITE(LuOut,80)
     &                                ((IT(I,J),J=1,3),I=1,3)
 80   FORMAT(3(1X,F10.5))
      RETURN
      END
