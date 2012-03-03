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
      SUBROUTINE GSCHMIDT(VEC,VORTH,NSIZE,NDIM,TMP,RESID,TOL)
C
C THIS PROJECTS OUT ALL PARTS OF AN INPUT VECTOR (VEC)
C WHICH LIE IN THE SPACE SPANNED BY THE ORTHOGONAL BASIS
C VORTH.
C
C   |v'> = |v> - SUM <i|v> |i>
C                 i
C
C WHERE THE |i> ARE NORMALIZED BASIS VECTORS FOR THE SPACE VORTH
C
C INPUT:
C       VEC : THE VECTOR WHICH IS TO BE ORTHOGONALIZED TO
C             THE EXISTING BASIS.  *IT IS ASSUMED THAT VEC
C             IS NORMALIZED ON INPUT)
C     VORTH : THE BASIS VECTORS FOR THE EXISTING ORTHOGONAL
C             BASIS
C     NSIZE : THE LENGTH OF THE BASIS VECTORS
C      NDIM : THE DIMENSION OF THE ORTHOGONAL SPACE
C       TMP : A SCRATCH VECTOR OF LENGTH NSIZE
C     RESID : THE NORM OF VEC, AFTER ORTHONALIZATION AND
C             BEFORE NORMALIZATION
C       TOL : TOLERANCE FOR RENORMALIZATION
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(NSIZE),VORTH(NSIZE,NDIM),TMP(NSIZE)
C
      DATA ONE /1.0/
C
      IF(NDIM.EQ.0)THEN
       RESID=XDNRM2(NSIZE,VEC,1)
       RETURN
      ENDIF
C
      CALL XDCOPY(NSIZE,VEC,1,TMP,1)
      DO 10 I=1,NDIM
       FACT=XDDOT(NSIZE,VORTH(1,I),1,VEC,1)
       CALL XDAXPY(NSIZE,-FACT,VORTH(1,I),1,TMP,1)
10    CONTINUE
      CALL XDCOPY(NSIZE,TMP,1,VEC,1)
C
C RENORMALIZE THE RESIDUAL
C
      RESID=XDNRM2(NSIZE,VEC,1)
      IF(RESID.GT.TOL)THEN
       X=ONE/RESID
       CALL XDSCAL(NSIZE,X,VEC,1)
      ENDIF 
      RETURN
      END
