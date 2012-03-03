      SUBROUTINE CUBEFORM (DENSACUB,DENSBCUB,NPTS,
     $     DENSTCUB,TMOACUB,TMOBCUB,NBASP,IDORO,IAORBD,
     $     IAORBO,TMOOUT,NBASPEFF,NORBOUTA,NORBOUTB,
     $     IORBOUT,IDENOUT)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION DENSACUB(NPTS+1)
      DIMENSION DENSBCUB(NPTS+1)
      DIMENSION DENSTCUB(NPTS+1)
      DIMENSION TMOACUB((NPTS+1)*NBASP)
      DIMENSION TMOBCUB((NPTS+1)*NBASP)
      DIMENSION TMOOUT((NPTS+1)*NBASPEFF)
C
C Calculate the Output density
      IF (IDORO.NE.1) THEN
         DO 10 I=1,(NPTS+1)
            IF (IAORBD.EQ.0) DENSTCUB(I)=DENSACUB(I)
            IF (IAORBD.EQ.1) DENSTCUB(I)=DENSBCUB(I)
            IF (IAORBD.EQ.2) DENSTCUB(I)=DENSACUB(I)+DENSBCUB(I)
C     WRITE(*,*) "output den #",I,"= ",DENSTCUB(I)
 10      CONTINUE
      END IF
C     
      NGRPS1=(NPTS+1)/6
      MODVAL1=MOD(NPTS+1,6)
      NGRPS2=(NPTS+1)*NBASPEFF/6
      MODVAL2=MOD((NPTS+1)*NBASPEFF,6)
C
C Write out the total density matrix
      IF (IDORO.NE.1) THEN
C
         DO 30 K=1,NGRPS1
            WRITE(IDENOUT,1000)
     $           (DENSTCUB((K-1)*6+JJ),JJ=1,6)
 1000       format(6(D13.5))
 30      CONTINUE
C     
         IF (MODVAL1.NE.0) THEN
            WRITE(IDENOUT,1001)
     $           (DENSTCUB(NGRPS1*6+JJ),JJ=1,MODVAL1)
 1001       format(6(D13.5))
         END IF
C     
      END IF
C     
C Construct the final orbital vector
      IF (IDORO.NE.0) THEN
C
         KK=0
         IF (IAORBO.EQ.0) THEN
           DO 50 K=1,((NPTS+1)*NBASP)
              TMOOUT(K)=TMOACUB(K)
 50        CONTINUE
        ELSE IF (IAORBO.EQ.1) THEN
           DO 51 K=1,((NPTS+1)*NBASP)
              TMOOUT(K)=TMOBCUB(K)
 51        CONTINUE
        ELSE
           DO 52 K=1,(NPTS+1)
              DO 53 KJ=1,NORBOUTA
                 KK=KK+1
                 TMOOUT(KK)=TMOACUB((K-1)*(NORBOUTA)+KJ)
 53           CONTINUE
C     
              DO 54 KJ=1,NORBOUTB
                 KK=KK+1
                 TMOOUT(KK)=TMOBCUB((K-1)*(NORBOUTB)+KJ)
 54           CONTINUE
 52        CONTINUE
        END IF
C
C     Write out the orbital coeffs.
        DO 60 K=1,NGRPS2
           WRITE(IORBOUT,1002)
     $          (TMOOUT((K-1)*6+JMO),JMO=1,6)
 1002      format(6(D13.5))
 60     CONTINUE
C     
        IF (MODVAL2.NE.0) THEN
           WRITE(IORBOUT,1003)
     $          (TMOOUT(NGRPS2*6+JMO),JMO=1,MODVAL2)
 1003      format(6(D13.5))
        END IF
C     
      END IF
C     WRITE(*,*)
C     WRITE(*,*)
C     
      RETURN
      END SUBROUTINE CUBEFORM

