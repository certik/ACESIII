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
      SUBROUTINE NORMAL(X,N)
C
C NORMALIZES VECTOR X
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)
      Q=0.D0
      DO 10 I=1,N
10    Q=Q+X(I)**2
      P=DSQRT(Q)
      IF(P.LT.1D-14)THEN
      WRITE(6,*)' null vector returned from NORMAL'
      RETURN
      ENDIF
      DO 11 I=1,N
11    X(I)=X(I)/P    
      RETURN
      END 
