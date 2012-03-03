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

      SUBROUTINE DBND(CARTCOORD,M,N,IBNDS,MAXREDUNCO,NRATMS,
     &                DERBMAT,IREDUNCO,TOTREDNCO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION LU
      INTEGER TOTREDNCO,A,B,M,N,IBNDS,NRATMS,IREDUNCO

      DIMENSION CARTCOORD(3*NRATMS),VECU(3),IREDUNCO(4,MAXREDUNCO),
c    &          DERBMAT(9*NRATMS*NRATMS*TOTREDNCO),
     &          DERBMAT(3*NRATMS,3*NRATMS*TOTREDNCO)

      DATA MONE/-1/

c     PRINT *,'*****ENTERING DBOND*****'
      LU = DIST(CARTCOORD(3*N-2),CARTCOORD(3*M-2))
      CALL VEC(CARTCOORD(3*M-2),CARTCOORD(3*N-2),VECU,1)
CSSS      print *,'BONDS'
CSSS      print *,'vecu',vecu(1),vecu(2),vecu(3)

      DO 20 IA=1,2
         A=IREDUNCO(IA,IBNDS)
         DO 30 IB=1,2
            B=IREDUNCO(IB,IBNDS)
            DO 40 I=1,3
            DO 50 J=1,3
               T1=(MONE**DELTA(A,B))*(VECU(I)*VECU(J)-DELTA(I,J))/LU
c               DERBMAT((3*B-3)+J,3*NRATMS*(IBNDS-1)+(3*A-3)+I)=T1
               DERBMAT((3*A-3)+I,3*NRATMS*(IBNDS-1)+(3*B-3)+J)=T1
50          CONTINUE
40          CONTINUE
30       CONTINUE
20    CONTINUE

CSSS       PRINT *,'DERIVATIVE OF THE BMAT FOR ATOMS',M,O,N
CSSS       CALL OUTPUT(DERBMAT,1,3*NRATMS,3*NRATMS*(IBNDS-1)+1,
CSSS     &              3*NRATMS*(IBNDS-1)+3*NRATMS,3*NRATMS,
CSSS     &              3*NRATMS*TOTREDNCO,1)

      RETURN
      END

