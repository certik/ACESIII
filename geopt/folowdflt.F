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

      SUBROUTINE FOLOWDFLT(HESMOD, DIAGHES, SCRATCH, VEC, TS, NRORMANR,
     &                     RFA, IVEC, IMODE, NCYCLE, NX, NOPT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL TS, NRORMANR, RFA

#include "io_units.par"

      DIMENSION HESMOD(NOPT, NOPT), DIAGHES(NOPT, NOPT),
     &          SCRATCH(NX*NX), VEC(NOPT)

      DATA IONE /1/

#include "jodaflags.com"

      CALL IGETREC(-1,'JOBARC','FIRSTCAL',IONE,ICYCLE)

      INEG = 0
      DO 1101 J=1, MAX(IVEC, 1)
         Z=100.D0
         DO 1100 I = 1, NOPT
            IF (SCRATCH(I).LT.Z) THEN
c              Assignment of eigenvector to be followed
               IF (ICYCLE.EQ.0.AND.TS) IMODE = I
               IMODE = I
               Z = SCRATCH(I)
            ENDIF
            IF (J.EQ.1.AND.HESMOD(I,I).LT.0.D0) INEG = INEG + 1
 1100    CONTINUE
         IF (ICYCLE.EQ.0.AND.TS) SCRATCH(IMODE) = 99.0D00
 1101 CONTINUE

C We follow nothing for minimas. Do this just to be on safe grounds.
      IF (NRORMANR) IMODE = 0

C On latter passes, determine overlap between Hessian eigenvectors
C and VEC (saved from the previous run). I am not sure whether this
C is the only strategy we should follow. See above for explanation.
C Ajith Perera 07/04

      IF (ICYCLE.GT.0.AND.TS) THEN
         ZOVLP = 0.d0
         DO 1833 I=1, NOPT
            SVH = xdot(NOPT,VEC,1,DIAGHES(1,I),1)
            IF (DABS(SVH).GT.ZOVLP) THEN
               ZOVLP = DABS(SVH)
               IMODE = I
            END IF
 1833    CONTINUE
         WRITE(LUOUT, 1351)IMODE, ZOVLP
 1351    FORMAT(T3,' Eigenvector ',i2,' has largest overlap with last ',
     &          'mode followed.',/,T3,' Magnitude of overlap is ',
     &           f8.5,'.')
      END IF
#ifdef _DEBUG_LVLM1

      WRITE(LUOUT, 2102) INEG
 2102 FORMAT(T3,' There are ',i2,' Negative Eigenvalues.')
#endif

      IF (INEG.GT.0.AND.NRORMANR) THEN
         IDIE = IFLAGS(50)
         IF (IDIE.EQ.0) then
            Write (LUOUT, 3131)
 3131       FORMAT(' @ANLYSHES, Negative eigenvalues in Hessian.')
            Call ERREX
         ElSE IF (IDIE.EQ.2) THEN
            WRITE(LUOUT,7141)
 7141       FORMAT(T3,' Negative eigenvalues are present in the ',
     $           'Hessian.'/t3,' Mode switched to rational function ',
     $           'approximation.')
            RFA = .TRUE.
            NRORMANR = .FALSE.
            IMODE = 0
         ELSE IF (IDIE.EQ.1) THEN
            DO 2411 I=1, NOPT
               HESMOD(I,I) = DABS(HESMOD(I,I))
 2411       CONTINUE
            WRITE(LUOUT,1241)
 1241       FORMAT(T3,' Negative eigenvalues present in Hessian have ',
     &         'been converted to their absolute values.')
         END IF
      END IF

c Save the eigenvector being followed for use in the next step
      DO 10 I = 1, NOPT
         VEC(I) = DIAGHES(I, IMODE)
 10   CONTINUE

C Now clean up the scratch array for future use.
      CALL IPUTREC(1,'JOBARC','FIRSTCAL',IONE,NCYCLE)
      CALL ZERO(SCRATCH, NX*NX)

      RETURN
      END

