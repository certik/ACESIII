
C CALCULATES THE VECTOR V BETWEEN CARTESIAN POINTS A AND B

      SUBROUTINE VEC(A,B,V,IX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),B(3),V(3)
      V(1)=B(1)-A(1)
      V(2)=B(2)-A(2)
      V(3)=B(3)-A(3)
      IF(IX.EQ.1)CALL NORMAL(V,3)
      RETURN
      END 
