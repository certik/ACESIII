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
