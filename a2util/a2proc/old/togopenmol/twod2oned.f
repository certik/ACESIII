      SUBROUTINE TWOD2ONED(PUTIN,NUMROWS,NUMCOLS,OUTPUT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
C     Takes a 2D matrix, reads down each column to
C     change it to a 1A array.
C     Ex: (1,1)->(2,1)->(3,1)->(1,2)->(2,2)->etc.
C           1      2      3      4      5
C
      DIMENSION PUTIN(NUMROWS,NUMCOLS)
      DIMENSION OUTPUT(NUMROWS*NUMCOLS)
C     
      DO 10 I=1,NUMROWS
         DO 20 II=(NUMCOLS-1),0,-1
            OUTPUT(I*NUMCOLS-II)=PUTIN(NUMCOLS-II,I)
            WRITE(*,*) "OUTPUT(",(I*NUMCOLS-II),")=PUTIN(",(NUMCOLS-II),
     $           ",",I,")=",PUTIN(NUMCOLS-II,I)
 20      CONTINUE
 10   CONTINUE
C     
      RETURN
      END SUBROUTINE TWOD2ONED
      
      
      
