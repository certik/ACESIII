      SUBROUTINE PRINTDVEC(DVECNAM,LEN)
C
C Designed to print out an integer vector
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION DVECNAM(LEN)
C     
      WRITE(*,*) "len=",LEN
      DO 20 I=1,LEN
         WRITE(*,*) DVECNAM(I)
 20   CONTINUE
C     
      END SUBROUTINE PRINTDVEC
      
      
