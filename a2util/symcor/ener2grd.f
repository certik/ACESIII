
C THIS ROUTINE CONSTRUCTS THE GRADIENT VECTOR FOR
C BLOCK NDIM FROM NUMERICAL DIFFERENTIATION OF THE ENERGY.

c INPUT
c integer NDIM
c double  ENERGY(2,NDIM)
c double  STPSIZ

c OUTPUT
c double GRD(NDIM) : symmetry coordinate gradient

      SUBROUTINE ENER2GRD(NDIM,ENERGY,GRD,STPSIZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ENERGY(2,NDIM), GRD(NDIM)

      LOGICAL PRINTQ

      COMMON /FLAGS/ IFLAGS(100)

      PRINTQ=(IFLAGS(1).GE.10)

c   o calculate gradient
      DTMP = 0.5d0/STPSIZ
      DO I=1,NDIM
         GRD(I)=(ENERGY(1,I)-ENERGY(2,I))*DTMP
      END DO

      IF (PRINTQ) THEN
         WRITE(6,1000)
1000     FORMAT(T3,'Internal coordinate gradient',/,
     &          T2,'Coord.',T16,'dE/dq')
         DO I=1,NDIM
            WRITE(6,'(2X,I5,5X,F20.10)')I,GRD(I)
         END DO
      END IF

      RETURN
      END

