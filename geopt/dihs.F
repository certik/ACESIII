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

      SUBROUTINE DDIH(CARTCOORD,N,P,O,M,IDIHS,TOTREDNCO,
     &                NRATMS,IREDUNCO,MAXREDUNCO,DERBMAT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION VECU,VECW,VECV,LU,LV,LW
      INTEGER A,B,I,J,TOTREDNCO,M,N,O,P,NRATMS
      DIMENSION CARTCOORD(3*NRATMS), IREDUNCO(4,MAXREDUNCO),
     &          VECU(3),VECW(3),VECV(3),PRODUW(3),PRODVW(3),
     &          DERBMAT(3*NRATMS,3*NRATMS*TOTREDNCO)
C
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0

CSSS      Write(6,"(4I5)") N, P, O, M
CSSS      Write(6,*)
      CALL VEC(CARTCOORD(3*O - 2),CARTCOORD(3*M - 2),VECV,1)
      CALL VEC(CARTCOORD(3*P - 2),CARTCOORD(3*O - 2),VECW,1)
      CALL VEC(CARTCOORD(3*P - 2),CARTCOORD(3*N - 2),VECU,1)
C
      ANGLUW = ANGLE(VECU, VECW, 3)
      ANGLVW = ANGLE(VECV, VECW, 3) 

CSSS       Write(6, "(2F10.5)") ANGLUW, ANGLVW 
       IF ((ANGLUW .EQ. 0.0D0 .OR. ANGLUW .GT. 175.0D0) .OR. 
     &     (ANGLVW .EQ. 0.0D0 .OR. ANGLWV .GT. 175.0D0)) RETURN
C
c      print *,'vecu',vecu(1),vecu(2),vecu(3)
c      print *,'vecv',vecv(1),vecv(2),vecv(3)
c      print *,'vecw',vecw(1),vecw(2),vecw(3)

      LU = DIST(CARTCOORD(3*O - 2), CARTCOORD(3*M - 2))
      LW = DIST(CARTCOORD(3*O - 2), CARTCOORD(3*P - 2))
      LV = DIST(CARTCOORD(3*P - 2), CARTCOORD(3*N - 2))

      CPHIU = DDOT(3,VECU,1,VECW,1)
      SPHIU2 =1.0D0-CPHIU*CPHIU
      CPHIV = DDOT(3,VECV,1,VECW,1)
      SPHIV2 =1.0D0-CPHIV*CPHIV

      CALL CROSS(VECU,VECW,PRODUW,0)
      CALL CROSS(VECV,VECW,PRODVW,0)

      D1 = (LU**2)*SPHIU2**2
      D2 = (LV**2)*SPHIV2**2
      D3 = 2.0D0*LU*LW*SPHIU2**2
      D4 = 2.0D0*LV*LW*SPHIV2**2
      D5 = 2.0D0*(LW**2)*SPHIU2**2
      D6 = 2.0D0*(LW**2)*SPHIV2**2
      D7 = LU*LW*SPHIU2
      D8 = LV*LW*SPHIV2

      DO 40 I=1,3
      DO 50 J=1,3
         T1I1 = (PRODUW(I)*(VECW(J)*CPHIU-VECU(J)))/D1
         T1I2 = (PRODUW(J)*(VECW(I)*CPHIU-VECU(I)))/D1
         T1   = T1I1+T1I2

         T2I1 = (PRODVW(I)*(VECW(J)*CPHIV-VECV(J)))/D2
         T2I2 = (PRODVW(J)*(VECW(I)*CPHIV-VECV(I)))/D2
         T2   = (T2I1+T2I2)

         T3I1 = PRODUW(I)*(VECW(J)-2.0D0*VECU(J)*
     &          CPHIU+VECW(J)*CPHIU**2)/D3
         T3I2 = PRODUW(J)*(VECW(I)-2.0D0*VECU(I)*
     &          CPHIU+VECW(I)*CPHIU**2)/D3
         T3   = (T3I1+T3I2)


         T4I1 = PRODVW(I)*(VECW(J)-2.0D0*VECV(J)*
     &          CPHIV+VECW(J)*CPHIV**2)/D4
         T4I2 = PRODVW(J)*(VECW(I)-2.0D0*VECV(I)*
     &          CPHIV+VECW(I)*CPHIV**2)/D4
         T4   = (T4I1+T4I2)

         T5I1 = (PRODUW(I)*(VECU(J)+(VECU(J)*CPHIU**2)-
     &          (3.0D0*VECW(J)*CPHIU)+(VECW(J)*CPHIU**3)))/D5
         T5I2 = (PRODUW(J)*(VECU(I)+(VECU(I)*CPHIU**2)-
     &          (3.0D0*VECW(I)*CPHIU)+(VECW(I)*CPHIU**3)))/D5
         T5   = (T5I1+T5I2)

         T6I1 = (PRODVW(I)*(VECV(J)+(VECV(J)*CPHIV**2)-
     &          (3.0D0*VECW(J)*CPHIV)+(VECW(J)*CPHIV**3)))/D6
         T6I2 = (PRODVW(J)*(VECV(I)+(VECV(I)*CPHIV**2)-
     &          (3.0D0*VECW(I)*CPHIV)+(VECW(I)*CPHIV**3)))/D6
         T6   = (T6I1+T6I2)

         K=0
         IF (I.NE.J) THEN
            IF      (I+J.EQ.3) THEN
               K = 3
            ELSE IF (I+J.EQ.4) THEN
               K = 2
            ELSE IF (I+J.EQ.5) THEN
               K = 1
            END IF
         END IF

         GAMMIJ = (J-I)*(-0.5D0)**(ABS(J-I))

         T7FACT=GAMMIJ
         T8FACT=GAMMIJ

         T7I = (VECW(K)*CPHIU-VECU(K))/D7
         T8I = (VECW(K)*CPHIV-VECV(K))/D8
         T7  = T7FACT*T7I
         T8  = T8FACT*T8I

C  DIAGONAL TERMS
C  11
       DERBMAT((O-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)=
     &   -t2-2.0D0*t4-t5+t6
C  22
       DERBMAT((P-1)*3+I,3*NRATMS*(IDIHS-1)+(P-1)*3+J)=
     &   t1-2.0D0*t3-t5+t6
C  33
       DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(N-1)*3+J)=
     &   t1
C  44
       DERBMAT((M-1)*3+I,3*NRATMS*(IDIHS-1)+(M-1)*3+J)=
     &   -t2

C  OFFDIAGONAL TERMS
C  21
       DERBMAT((P-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)=
     &    t3+t4+t5-t6+t7-t8

      DERBMAT((O-1)*3+J,3*NRATMS*(IDIHS-1)+(P-1)*3+I)=
     & DERBMAT((P-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)

C  31
       DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)=
     &    -t3-t7

       DERBMAT((O-1)*3+J,3*NRATMS*(IDIHS-1)+(N-1)*3+I)=
     & DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)

C  41
       DERBMAT((M-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)=
     &    t2+t4+t8

       DERBMAT((O-1)*3+J,3*NRATMS*(IDIHS-1)+(M-1)*3+I)=
     & DERBMAT((M-1)*3+I,3*NRATMS*(IDIHS-1)+(O-1)*3+J)

C  32
       DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(P-1)*3+J)=
     &    -t1+t3+t7

       DERBMAT((P-1)*3+J,3*NRATMS*(IDIHS-1)+(N-1)*3+I)=
     & DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(P-1)*3+J)

C  42
       DERBMAT((M-1)*3+I,3*NRATMS*(IDIHS-1)+(P-1)*3+J)=
     &    -t4-t8

       DERBMAT((P-1)*3+J,3*NRATMS*(IDIHS-1)+(M-1)*3+I)=
     & DERBMAT((M-1)*3+I,3*NRATMS*(IDIHS-1)+(P-1)*3+J)

C  34
       DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(M-1)*3+J)=
     &    0.0D0

c      DERBMAT((M-1)*3+J,3*NRATMS*(IDIHS-1)+(N-1)*3+I)=
c    & DERBMAT((N-1)*3+I,3*NRATMS*(IDIHS-1)+(M-1)*3+J)

c      PRINT *,'I,J,K',I,J,K
c      print *,'33',t1 - 2.0D0*t3 - t5 + t6
c      PRINT *,'23',-t1 + t3 + GAMMIJ*t7I
c      print *,'34',t3 + t4 + t5 - t6
c    &                      + GAMMIJ*(t7I-T8I)
c      print *,'44',-t2 - 2.0D0*t4 + t6 - t5
c      print *,'24',-t3 - GAMMIJ*t7I
c      print *,'53',-t4 - GAMMIJ*T8I
c      print *,'54',t2 + t4 + GAMMIJ*T8I
c      print *,'55',-t2
c      PRINT *,'22',t1

50    CONTINUE
40    CONTINUE

CSSS      PRINT *,'DERIVATIVE OF THE BMAT FOR ATOMS',N,P,O,M
CSSS      CALL OUTPUT(DERBMAT,1,3*NRATMS,3*NRATMS*(IDIHS-1)+1,
CSSS     &              3*NRATMS*(IDIHS-1)+3*NRATMS,3*NRATMS,
CSSS     &              3*NRATMS*TOTREDNCO,1)

      RETURN
      END

