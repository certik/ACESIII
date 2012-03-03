C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.

      SUBROUTINE DEANG(CARTCOORD,M,O,N,BMATRX,NRATMS,TOTREDNCO,IREDUNCO,
     &                MAXREDUNCO,IANGS,DERBMAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DOUBLE PRECISION LU,LV
      INTEGER A,B,M,O,N,NRATMS,TOTREDNCO,IREDUNCO,MAXREDUNCO,IANGS
      DIMENSION CARTCOORD(3*NRATMS),VECU(3),VECV(3),
     &          IREDUNCO(4,MAXREDUNCO),BMATRX(TOTREDNCO,3*NRATMS),
     &          DERBMAT(3*NRATMS,3*NRATMS*TOTREDNCO)

      CALL VEC(CARTCOORD(3*O-2),CARTCOORD(3*M-2),VECU,1)
      CALL VEC(CARTCOORD(3*O-2),CARTCOORD(3*N-2),VECV,1)

      LU = DIST(CARTCOORD(3*O-2),CARTCOORD(3*M-2))
      LV = DIST(CARTCOORD(3*O-2),CARTCOORD(3*N-2))
      COSQ= DDOT(3,VECU,1,VECV,1)
      SINQ= DSQRT(1-(COSQ**2))
      D1 = (LU**2)*SINQ
      D2 = (LV**2)*SINQ
      D3 = (LU*LV)*SINQ

CSSS      print *,'*****Entering DEANG*****'
CSSS      print *,'vecu',vecu(1),vecu(2),vecu(3)
CSSS      print *,'vecv',vecv(1),vecv(2),vecv(3)
CSSS      PRINT *,'MON',M,O,N
CSSS      Print *, "SINQ,COSQ", SINQ, COSQ

      DO 10 I=1,3
         BM1=BMATRX(IANGS,(3*O-3)+I)
         BM2=BMATRX(IANGS,(3*M-3)+I)
         BM3=BMATRX(IANGS,(3*N-3)+I)
      DO 20 J=1,3
         BM4=BMATRX(IANGS,(3*O-3)+J)
         BM5=BMATRX(IANGS,(3*M-3)+J)
         BM6=BMATRX(IANGS,(3*N-3)+J)
C MM TERM
         BMMM=(COSQ/SINQ)*BM5*BM2
C OO TERM
         BMOO=(COSQ/SINQ)*BM4*BM1
C NN TERM
         BMNN=(COSQ/SINQ)*BM6*BM3
C MO TERM
         BMMO=(COSQ/SINQ)*BM4*BM2
C NO TERM
         BMNO=(COSQ/SINQ)*BM4*BM3
C MN TERM
         BMMN=(COSQ/SINQ)*BM6*BM2

         T1= (VECU(I)*VECV(J)+VECU(J)*VECV(I)-3.0D0*VECU(I)
     &       *VECU(J)*COSQ+DELTA(I,J)*COSQ)/D1
         T2= (VECV(I)*VECU(J)+VECV(J)*VECU(I)-3.0D0*VECV(I)
     &       *VECV(J)*COSQ+DELTA(I,J)*COSQ)/D2
         T3= (VECU(I)*VECU(J)+VECV(J)*VECV(I)-VECU(I)*VECV(J)*COSQ-
     &       DELTA(I,J))/D3
         T4= (VECV(I)*VECV(J)+VECU(J)*VECU(I)-VECV(I)*VECU(J)*COSQ-
     &       DELTA(I,J))/D3

C Diagonal Terms
C  11
         DERBMAT((O-1)*3+I,3*NRATMS*(IANGS-1)+(O-1)*3+J)=
     &          t1+t2+t3+t4-BMOO
C  22
         DERBMAT((M-1)*3+I,3*NRATMS*(IANGS-1)+(M-1)*3+J)=
     &          t1-BMMM
C  33
         DERBMAT((N-1)*3+I,3*NRATMS*(IANGS-1)+(N-1)*3+J)=
     &          t2-BMNN

C  Offdiagonal Term
C  23
         DERBMAT((M-1)*3+I,3*NRATMS*(IANGS-1)+(N-1)*3+J)=
     &          t3-BMMN

         DERBMAT((N-1)*3+J,3*NRATMS*(IANGS-1)+(M-1)*3+I)=
     &   DERBMAT((M-1)*3+I,3*NRATMS*(IANGS-1)+(N-1)*3+J)
C  21
         DERBMAT((M-1)*3+I,3*NRATMS*(IANGS-1)+(O-1)*3+J)=
     &           -t1-t3-BMMO
         DERBMAT((O-1)*3+J,3*NRATMS*(IANGS-1)+(M-1)*3+I)=
     &   DERBMAT((M-1)*3+I,3*NRATMS*(IANGS-1)+(O-1)*3+J)
C  31
         DERBMAT((N-1)*3+I,3*NRATMS*(IANGS-1)+(O-1)*3+J)=
     &           -t2-t4-BMNO
         DERBMAT((O-1)*3+J,3*NRATMS*(IANGS-1)+(N-1)*3+I)=
     &   DERBMAT((N-1)*3+I,3*NRATMS*(IANGS-1)+(O-1)*3+J)

  20  CONTINUE
  10  CONTINUE

CSSS       print *,'DERIVATIVE OF THE BMAT FOR ATOMS',M,O,N
CSSS       CALL OUTPUT(DERBMAT,1,3*NRATMS,3*NRATMS*(IANGS-1)+1,
CSSS     &              3*NRATMS*(IANGS-1)+3*NRATMS,3*NRATMS,
CSSS     &              3*NRATMS*TOTREDNCO,1)

      RETURN
      END

