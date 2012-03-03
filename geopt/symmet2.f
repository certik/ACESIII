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

C THIS ROUTINE SYMMETRIZES A GIVEN MATRIX A
C
C   A(PQ) = 1/2 ( A(PQ)+A(QP))
C
C WHERE A IS A SYMMETRY PACKED MATRIX AND
C NUM THE CORRESPONDING POPULATION VECTOR
C
C THE SYMMETRIZATION IS HERE COMPLETELY DONE IN PLACE
C
C CODED AUGUST/90 JG

      SUBROUTINE SYMMET2(A,NUM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NUM,NUM)

      if (num.lt.2) return

C SYMMETRIZE FIRST HALF OF A
      DO I = 2, NUM
         DO J = 1, I-1
            A(I,J) = 0.5d0*(A(I,J)+A(J,I))
         END DO
      END DO

C FILL SECOND HALF OF A
      DO I = 2, NUM
         DO J = 1, I-1
            A(J,I) = A(I,J)
         END DO
      END DO

      RETURN
      END
