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
      SUBROUTINE CURVE(SCRATCH, EHESS, CURVTRE, NDIM, NX)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION SCRATCH(NX*NX), EHESS(NDIM, NDIM)
C
      DATA SUM1 /0.0D0/, SUM2 /0.0D0/, SUM3 /0.0D0/
      
      DO 10 I = 1, NDIM
C
         SUM1 = SUM1 + EHESS(I, I)**4*SCRATCH(I + 2*NDIM)**2
         SUM2 = SUM2 + EHESS(I, I)**2*SCRATCH(I + 2*NDIM)**2
         SUM3 = SUM3 + EHESS(I, I)**3*SCRATCH(I + 2*NDIM)**2
C
 10   CONTINUE
      
      CURVTRE = DSQRT(DABS(SUM1*SUM2 - SUM3**2))/DSQRT(DABS((SUM2**3)))
C
      RETURN
      END
