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

C SUBROUTINE COMPUTES OUTER PRODUCT OF VECTOR V WITH ITSELF.  IF
C ICOMP IS SET TO ONE, 1-|V><V| IS RETURNED; OTHERWISE IT IS JUST
C |V><V|.

      SUBROUTINE OUTERP(VV,V,NR,NC,ICOMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(NR,NC),VV(NR,NR)
      DO I = 1, NR
         DO J = 1, NR
            T = 0.D0
            DO K = 1, NC
               T = T + V(I,K) * V(J,K)
            END DO
            VV(I,J) = T
         END DO
      END DO
      IF (ICOMP.EQ.0) RETURN
      Z = -1.D0
      CALL xscal(NR*NR,Z,VV,1)
      DO I = 1, NR
         VV(I,I) = 1.D0 + VV(I,I)
      END DO
      RETURN
      END

