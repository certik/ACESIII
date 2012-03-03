      SUBROUTINE PRINTMAT(VECNAM,NUMROWS,NUMCOLS)
C     
C     Designed to print out a double-precision
c     matrix from its vector form
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION VECNAM(NUMROWS*NUMCOLS)
C     
C     DO 20 I=1,(NUMROWS*NUMCOLS)
C        WRITE(*,*) VECNAM(I)
C     20 CONTINUE
C     
      NGRPS=NUMCOLS/5 
      MODVAL=MOD(NUMCOLS,5)
      DO 5 J=0,(NGRPS-1)
         DO 10 I=1,NUMROWS
            WRITE(*,"(5 4D12.4)") (VECNAM(J*
     $           NUMROWS*5+NUMROWS*II+I),II=0,4)
 10      CONTINUE
         WRITE(*,*)
         WRITE(*,*)
 5    CONTINUE
      IF (MODVAL.NE.0) THEN
         DO 20 I=1,NUMROWS
            WRITE(*,"(5 4D12.4)") (VECNAM(NGRPS*
     $           NUMROWS*5+NUMROWS*II+I),II=0,
     $           MODVAL-1)
 20      CONTINUE
         WRITE(*,*)
         WRITE(*,*)
      END IF
C     
      END SUBROUTINE PRINTMAT
      
