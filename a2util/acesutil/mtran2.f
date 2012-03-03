
C IN-PLACE TRANSPOSITION OF AN N x N MATRIX.

      SUBROUTINE MTRAN2(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N)
      if (n.lt.1) return
      DO I = 1, N
         DO J = 1, I-1
            Z      = A(I,J)
            A(I,J) = A(J,I)
            A(J,I) = Z
         END DO 
      END DO
      RETURN
      END
