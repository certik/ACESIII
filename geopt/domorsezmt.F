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
      SUBROUTINE DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C     Maximum number of atoms currently allowed
#include "mxatms.par"
#include "cbchar.com"
#include "coord.com"
C
      DIMENSION SCRATCH(NX*NX)
C
      DO 1452 I = 1, NOPT
C
         IIN = 0
C     
         DO 1453 IJ = 1, NATOMS
C
C NOPTI(I) gives "SQUASHED" position - need unsquashed for connectivity
C
            IK = NOPTI(I) + 6 - 3/NOPTI(I)
            IF(IJ*3 + 1 .EQ. IK) IIN = IK
C
 1453    CONTINUE
C
         IF (IIN .NE. 0) THEN
            IAT = 1 + (IIN - 1)/3
            JAT = NCON(IIN)
            NA = IATNUM(IAT)
            NB = IATNUM(JAT)
C
            CALL MORSEA(NA, NB, A)
C     
            IF(A.EQ.0.D0)THEN
               WRITE(LUOUT,1474)ZSYM(IAT),ZSYM(JAT)
 1474          FORMAT(T3,'@EFOL-I, No Morse constant available',
     &                'for ',A5, '-',A5,'.  Default used.')
            ENDIF
C
C Scale factor calculated from internal coordinates (not
C symmetry coordinate increment)
C
            Z = SCRATCH(NOPT + I)/DSQRT(DFLOAT(NEQ(NOPTI(I)) + 1))
C
            IF (DABS(Z) .GT. 0.75D0) THEN
               WRITE(LUOUT, 7208) VARNAM(IK)
 7208          FORMAT(T3,' NR step for ',A5,' too large.  MANR ',
     &                'scaling not done.')
            ELSE
               CORR = 1.D0 + 1.5D0*A*Z+2.3333333D0*(A*Z)**2
               WRITE(LUOUT, 7201)CORR,VARNAM(IK)
 7201          FORMAT(T3,' MANR scale factor for NR step is ',
     &                F5.3,' for ', A5,'.')
               SCRATCH(NOPT + I) = SCRATCH(NOPT + I)*CORR
            ENDIF
C     
         ENDIF
C
 1452 CONTINUE
C
      RETURN
      END


