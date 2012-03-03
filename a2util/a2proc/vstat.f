      SUBROUTINE VSTAT(V,ZQ,LENGTH)
C
C RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
C     ZQ(1)  Largest absolute magnitude
C     ZQ(2)  Smallest absolute magnitude
C     ZQ(3)  Largest value
C     ZQ(4)  Smallest value
C     ZQ(5)  2-norm
C     ZQ(6)  Dynamic range of the vector (abs. min. - abs. max.)
C   

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LENGTH),ZQ(6)
      U=0.D0
      CALL ZERO(ZQ,6)
      ZQ(2)=DABS(V(1))
      ZQ(4)=V(1)
      DO 20 I=1,LENGTH
      ZQ(1)=MAX(ZQ(1),DABS(V(I)))
      ZQ(2)=MIN(ZQ(2),DABS(V(I)))
      ZQ(3)=MAX(ZQ(3),V(I))
      ZQ(4)=MIN(ZQ(4),V(I))
20    U=U+V(I)*V(I)
      If (Length .ne. 0) ZQ(5)=DSQRT(U/LENGTH)
      ZQ(6)=ZQ(2)-ZQ(1)
      RETURN
      END
