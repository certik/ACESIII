      SUBROUTINE INITP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,EOR,DSTRT
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
      PARAMETER (LENK  = 8*MXCBF + 9 + 7 + 8 + 3 + 7*28 + 64 + 8*8*2
     &                 + 8*2 + 100 + MXTNSH*8) 
C-----------------------------------------------------------------------
      COMMON /INDX/ PC(512),K(LENK)
c     COMMON /INDX/ PC(512),K(5011)
c     COMMON /INDX/ PC(512),K(512)
C     INITIALIZE THE POPULATION COUNT FACTOR
      K(1)=0
      N=1
      DO 1 I=1,9
        DO 2 J=1,N
    2   K(J+N) = K(J)+1
    1 N = N + N
      DO 3 I=1,512
    3 PC(I) = DFLOAT(1 - 2*AND(K(I),1))
      RETURN
      END
