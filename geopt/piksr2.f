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

C SIMPLE SORTER FROM NUMERICAL RECIPES.

      SUBROUTINE PIKSR2(N,ARR,NLIST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARR(N),NLIST(N)
      if (n.lt.2) return
      DO J = 2, N
         A = ARR(J)
         NLST = NLIST(J)
         DO I = J-1, 1, -1
            IF (ARR(I).LE.A) GOTO 10
            ARR(I+1) = ARR(I)
            NLIST(I+1) = NLIST(I)
         END DO
         I = 0
   10    ARR(I+1) = A
         NLIST(I+1) = NLST
      END DO
      RETURN
      END
