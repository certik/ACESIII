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
      SUBROUTINE ALONGX(Q1,Q2,RM,IERR)
C
C THIS SUBROUTINE RETURNS THE ROTATION MATRIX NEEDED TO ROTATE
C THE XY-PROJECTION OF THE BISECTOR OF THE TWO 3-VECTORS (Q1 AND
C Q2) ALONG THE X AXIS.  IERR=-1 MEANS THAT THE TWO VEC
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Q1(3),Q2(3),RM(9)
      DATA ONEM  /-1.0/
      DATA ONE   / 1.0/
      DATA ZILCH / 0.0/
      DATA TOL  /1.D-10/
C
      IERR=0
C
C CALCULATE DIFFERENCE VECTOR
C
      CALL VADD(RM,Q1,Q2,3,ONEM) 
      TMP=RM(1)
      RM(1)=RM(2)
      RM(2)=-TMP
      RM(3)=ZILCH
      X=XDNRM2(3,RM,1)
      IF(X.LT.TOL)THEN
       IERR=-1
       RETURN
      ELSE
       IERR=0
       FACT=ONE/X
       CALL DSCAL(3,FACT,RM,1)
      ENDIF
C
      ANGLE=SIGN(DACOSX(RM(1),1.D-10),RM(2))
      CALL ROTM(3,ANGLE,0,RM)
C
      RETURN
      END
