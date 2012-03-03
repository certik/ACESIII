
C IN-PLACE TRANSPOSITION OF A SUBMATRIX OF Z, WHERE Z IS
C AN LDIM x NDIM MATRIX.  THE UPPER NDIM x NDIM TRIANGLE
C IS TRANSPOSED IN THIS ROUTINE

      SUBROUTINE MTRAN3(Z,LDIM,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(LDIM,NDIM)
      if ((ldim.lt.1).or.(ndim.lt.1)) return
      DO J = 1, NDIM
         DO I = 1, J-1
            XTMP   = Z(I,J)
            Z(I,J) = Z(J,I)
            Z(J,I) = XTMP
         END DO
      END DO
      RETURN
      END
