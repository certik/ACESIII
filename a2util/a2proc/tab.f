C
      SUBROUTINE TAB(NOUT,A,N,M,NNN,MMM)
C
C  modified by RPM for an 80-character wide display.  un-comment the 
C  appropriate set of parameters and format statements to choose how
C  many columns across the page.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION A(NNN,MMM)
       integer wid
C WID 120
C      parameter(wid = 10)
C    1 FORMAT(1X,I3,3X,10F12.6)
C   12 FORMAT(4X,10(8X,I4))
C wid 80
C       parameter(wid=5)
C    1 FORMAT(1X,I3,3X,5e14.6)
C12    FORMAT(4X, 5(8x,i4))
C      parameter (wid = 4)
C   1 FORMAT(1X,I3,1X,4e18.10)
C12    FORMAT(4X, 4(8x,i4))
       parameter(wid=6)
 1     FORMAT(1X,I3,2X,6F12.9)
 12    FORMAT(4X, 6(8x,i4))
C
   11 FORMAT(/)
      MM=M/wid
      IF(MM.EQ.0) GO TO 6
          DO 4 II=1,MM
             JP=(II-1)*wid+1
             JK=II*wid
             WRITE(nout,11)
             WRITE(nout,12)(I,I=JP,JK)
             DO 4 I=1,N
                WRITE(nout,1)I,(A(I,J),J=JP,JK)
    4     CONTINUE
    6 CONTINUE
      MA=MM*wid+1
      IF(MA .GT. M) RETURN
      WRITE(nout,11)
      WRITE(nout,12)(I,I=MA,M)
      DO 5 I=1,N
          WRITE(nout,1) I,(A(I,J),J=MA,M)
    5 CONTINUE
      RETURN
      END
