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

      SUBROUTINE COMPKMAT(DERBMAT,NRATMS,TOTREDNCO,FI,HC,
     &                    DIFTEMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER TOTREDNCO
      DOUBLE PRECISION KMAT

#include "mxatms.par"

C KMAT is used as temporary storge while building DIFTEMP and it is of
C size 9*MXATMS*MXATMS (so it needs to be managed dynamically).

      DIMENSION DERBMAT(3*NRATMS,3*NRATMS*TOTREDNCO),
     &          FI(TOTREDNCO),HC(3*NRATMS,3*NRATMS),
     &          KMAT(3*NRATMS,3*NRATMS),
     &          DIFTEMP(3*NRATMS,3*NRATMS)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      CALL ZERO(KMAT,9*NRATMS*NRATMS)

      LENDBMAT = 9*NRATMS*NRATMS*TOTREDNCO
      CALL DGETREC(20,'JOBARC','DERBMAT',LENDBMAT,DERBMAT)

      DO I=1, TOTREDNCO
         IND=3*NRATMS*(I-1)+1
C
#ifdef _DEBUG_LVL0
         Write(6, "(a,4F10.5)") 'The gradients Compkmat: ',
     &        (FI(J), J=1, TOTREDNCO)
#endif
         CALL DSCAL(9*NRATMS*NRATMS,FI(I),DERBMAT(1,IND),1)
      END DO

      DO I=1,3*NRATMS
         DO J=1,3*NRATMS
            DO K=1,TOTREDNCO
               JTMP = (K-1)*3*NRATMS+J
               KMAT(I,J) = KMAT(I,J)+ DERBMAT(I,JTMP)
            END DO
         END DO
      END DO
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") 'The  Cartesian Hess'
      CALL OUTPUT(HC,1,3*NRATMS,1,3*NRATMS,3*NRATMS,
     &            3*NRATMS,1)
#endif
C
      DO I=1,3*NRATMS
         DO J=1,3*NRATMS
            DIFTEMP(I,J)=HC(I,J)-KMAT(I,J)
         END DO
      END DO
C
#ifdef _DEBUG_LVL0
      Write(6,*) 'In COMPKMAT: The (HC-K) Matrix'
      CALL OUTPUT(DIFTEMP,1,3*NRATMS,1,3*NRATMS,3*NRATMS,
     &            3*NRATMS,1)
#endif

      RETURN
      END

