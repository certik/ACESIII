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
      SUBROUTINE CONVQ(A,FC,FI,IX,ISYM)
C
C CONVERTS FORCES IN CARTESIAN COORDINATES TO INTERNAL COORDINATES (IX=0
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL XYZIN, NWFINDIF
#include "mxatms.par"
#include "flags.h"
#include "jodaflags.com"
#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
          DIMENSION A(3*NATOMS,NXM6),FC(NX),FI(NXM6)
C
#include "io_units.par"
C
#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,*) "@-CONVQF, The incoming Cartesian gradients"
         WRITE(6,"(6F13.7)") (FC(IJL), IJL = 1, NX)
         Write(6,*) "@-CONVQF, The A matrix"
         CALL OUTPUT(A, 1,3*NATOMS, 1, NXM6, 3*NATOMS, NXM6, 1)
#endif
      IF(IX.EQ.0)THEN
         DO 10 I=1,NXM6
            Z=0.D0
            DO 20 J=1,NX
               Z=FC(J)*A(J,I)+Z
 20         CONTINUE
            FI(I)=Z
 10      CONTINUE

#ifdef _DEBUG_LVL0
         Write(6,*) "@-CONVQF,  The outgoing internal gradients"
         WRITE(6,*) (FI(IJL), IJL = 1, NXM6)
#endif
C
C SYMMETRIZE THE GRADIENT IF NECESSARY
C
         IF (ISYM.EQ.1)THEN
            IF (XYZIN) THEN
               IF (IFLAGS(h_IFLAGS_symmetry).EQ. 2) THEN
                  CALL IGETREC(20, 'JOBARC', 'PLSMINSP', NXM6, NCON)
               END IF
               DO I = 1, NUNIQUE
                  IP = IUNIQUE(I)
cSSS                  WRITE(6,*) "The IP =", IP, NEQ(IP),FI(IP) 
                  Z  = FI(IP)
C
                  DO J = 1, NEQ(IP)
                     IF (NCON(IEQUIV(I, J)) .EQ. 1) THEN
                        SIGN = -1.0D0
                     ELSE
                        SIGN =  1.0D0
                     ENDIF
cSSS                     WRITE(6,*) FI(IEQUIV(I,J)), IEQUIV(I,J),SIGN
                     Z = Z + SIGN*FI(IEQUIV(I, J))
                  ENDDO
C
                  FIAVE=Z/(NEQ(IP)+1)
cSSS                  WRITE(6,*) FIAVE
                  DIFF=DABS(FI(IP)-FIAVE)
                  IF(DIFF.GT.1D-7)WRITE(LuOut,72)DIFF,IP
C
 72               FORMAT(T3,'@CONVQ-W, Nominally equivalent',
     $            'internal coordinate gradients vary by at least ',
     $            e10.5,'.'/
     $            T3,' Problem with gradients equivalent to ',
     $            'parameter [',I3,'].')
                  FI(IP) = FIAVE
C
                  DO J = 1, NEQ(IP)
                     IF (NCON(IEQUIV(I, J)) .EQ. 1) THEN
                        SIGN = -1.0D0
                     ELSE
                        SIGN = 1.0D0
                     ENDIF
                     FI(IEQUIV(I, J))= SIGN*FIAVE
cSSS                     WRITE(6,*) "The IP =", IEQUIV(I,J),FI(IEQUIV(I, J))
                  ENDDO
               ENDDO
               CALL FILTER(FI, NXM6, 1.0D-6) 
               RETURN
            ELSE
               DO 100 I=1,NUNIQUE
                  IP=IUNIQUE(I)
                  Z=FI(IP)
                  DO 102 J=1,NEQ(IP)
 102                 Z=Z+FI(IEQUIV(I,J))
                     FIAVE=Z/(NEQ(IP)+1)
                     DIFF=DABS(FI(IP)-FIAVE)
                     IF(DIFF.GT.1D-7)WRITE(LuOut,71)DIFF,IP
 71                 FORMAT(T3,'@CONVQ-W, Nominally equivalent',
     $             'internal coordinate gradients vary by at least ',
     $              e10.5,'.'/
     $              T3,' Problem with gradients equivalent to ',
     $              'parameter [',I3,'].')
                    FI(IP)=FIAVE
C
                    DO J = 1, NEQ(IP)
                       FI(IEQUIV(I,J))= FIAVE
                    ENDDO
C
 100             CONTINUE
            ENDIF
            RETURN
         ENDIF
      ELSE
         DO 30 I=1,NX
            Z=0.D0
            DO 40 J=1,NXM6
 40            Z=A(I,J)*FI(J)+Z
            FC(I)=Z
 30      CONTINUE
         RETURN
      ENDIF
      RETURN
      END
