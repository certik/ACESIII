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

C THIS ROUTINE SETS ALL VALUES IN A VECTOR A TO ZERO IF
C THEY ARE BELOW A SPECIFIED TOLERANCE (TOL).

      SUBROUTINE FILTER(A,LENGTH,TOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LENGTH)
      if (length.lt.1) return
      DO I = 1, LENGTH
         X = ABS(A(I))
         IF (X.LT.TOL) A(I) = 0.0
      END DO
      RETURN
      END
