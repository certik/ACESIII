      SUBROUTINE LAME(NUMSHL,IXVEC,IYVEC,IZVEC,IORD)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IXVEC(IORD),IYVEC(IORD),IZVEC(IORD)
C
      IF (NUMSHL.EQ.1) THEN
         IXVEC(1)=0
         IYVEC(1)=0
         IZVEC(1)=0
      ELSEIF (NUMSHL.EQ.2) THEN
         IXVEC(1)=1
         IYVEC(1)=0
         IZVEC(1)=0
         IXVEC(2)=0
         IYVEC(2)=1
         IZVEC(2)=0
         IXVEC(3)=0
         IYVEC(3)=0
         IZVEC(3)=1
      ELSEIF (NUMSHL.EQ.3) THEN
         IXVEC(1)=2
         IYVEC(1)=0
         IZVEC(1)=0
         IXVEC(2)=1
         IYVEC(2)=1
         IZVEC(2)=0
         IXVEC(3)=1
         IYVEC(3)=0
         IZVEC(3)=1
         IXVEC(4)=0
         IYVEC(4)=2
         IZVEC(4)=0
         IXVEC(5)=0
         IYVEC(5)=1
         IZVEC(5)=1
         IXVEC(6)=0
         IYVEC(6)=0
         IZVEC(6)=2
      ELSEIF (NUMSHL.EQ.4) THEN
         IXVEC(1)=3
         IYVEC(1)=0
         IZVEC(1)=0
         IXVEC(2)=2
         IYVEC(2)=1
         IZVEC(2)=0
         IXVEC(3)=2
         IYVEC(3)=0
         IZVEC(3)=1
         IXVEC(4)=1
         IYVEC(4)=2
         IZVEC(4)=0
         IXVEC(5)=1
         IYVEC(5)=1
         IZVEC(5)=1
         IXVEC(6)=1
         IYVEC(6)=0
         IZVEC(6)=2
         IXVEC(7)=0
         IYVEC(7)=3
         IZVEC(7)=0
         IXVEC(8)=0
         IYVEC(8)=2
         IZVEC(8)=1
         IXVEC(9)=0
         IYVEC(9)=1
         IZVEC(9)=2
         IXVEC(1)=0
         IYVEC(1)=0
         IZVEC(1)=3
      ELSEIF (NUMSHL.EQ.5) THEN
         IXVEC(1)=4
         IYVEC(1)=0
         IZVEC(1)=0
         IXVEC(2)=3
         IYVEC(2)=1
         IZVEC(2)=0
         IXVEC(3)=3
         IYVEC(3)=0
         IZVEC(3)=1
         IXVEC(4)=2
         IYVEC(4)=2
         IZVEC(4)=0
         IXVEC(5)=2
         IYVEC(5)=1
         IZVEC(5)=1
         IXVEC(6)=2
         IYVEC(6)=0
         IZVEC(6)=2
         IXVEC(7)=1
         IYVEC(7)=3
         IZVEC(7)=0
         IXVEC(8)=1
         IYVEC(8)=2
         IZVEC(8)=1
         IXVEC(9)=1
         IYVEC(9)=1
         IZVEC(9)=2
         IXVEC(1)=1
         IYVEC(1)=0
         IZVEC(1)=3
         IXVEC(11)=0
         IYVEC(11)=4
         IZVEC(11)=0
         IXVEC(12)=0
         IYVEC(12)=3
         IZVEC(12)=1
         IXVEC(13)=0
         IYVEC(13)=2
         IZVEC(13)=2
         IXVEC(14)=0
         IYVEC(14)=1
         IZVEC(14)=3
         IXVEC(15)=0
         IYVEC(15)=0
         IZVEC(15)=4
      END IF
C
      RETURN
      END SUBROUTINE LAME

