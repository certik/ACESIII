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
      SUBROUTINE DOMORSEXYZ_RIC(SCRATCH, NOPT, NX, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER TOTREDNCO, TOTNOFBND
      LOGICAL BONDED
C     Maximum number of atoms currently allowed
#include "mxatms.par"
C
      DIMENSION IBNDTO(NATOMS*NATOMS), IREDUNCO(4, MAXREDUNCO) 
      DIMENSION SCRATCH(NX*NX)
      LOGICAL I_UNIQUE(MAXREDUNCO)
C
#include "cbchar.com"
#include "coord.com"
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C 
      CALL IGETREC(20, 'JOBARC', 'IBONDTO ', NATOMS*NATOMS, IBNDTO)
      CALL IGETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
      CALL IGETREC(20, 'JOBARC', 'UNIQUEDF', TOTREDNCO, I_UNIQUE)
      CALL IGETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO, IREDUNCO)
C
      DO IOPT = 1,  NOPT
C
         IATM = IREDUNCO(1, IOPT)
         JATM = IREDUNCO(2, IOPT)
         KATM = IREDUNCO(3, IOPT)
C
         IF (I_UNIQUE(IOPT) .AND. KATM .EQ. 0) THEN
C
            NA = IATNUM(IATM)
            NB = IATNUM(JATM)
C
            CALL MORSEA(NA, NB, A)
C     
            IF(A.EQ.0.D0)THEN
               WRITE(LUOUT,1474)ZSYM(IATM),ZSYM(IATM)
 1474          FORMAT(T3,'@EFOL-I, No Morse constant available',
     &                'for ',A5, '-',A5,'.  Default used.')
            ENDIF
C
C Scale factor calculated from internal coordinates (not
C symmetry coordinate increment)
C
            Z = SCRATCH(NOPT + IOPT)/
     &          DSQRT(DFLOAT(NEQ(NOPTI(IOPT)) + 1))
C
            IF (DABS(Z) .GT. 0.75D0) THEN
               WRITE(LUOUT, 7208) VARNAM(IOPT)
 7208          FORMAT(T3,' NR step for ',A5,' too large.  MANR ',
     &                'scaling not done.')
            ELSE
               CORR = 1.D0 + 1.5D0*A*Z+2.3333333D0*(A*Z)**2
               WRITE(LUOUT, 7201)CORR,VARNAM(IOPT)
 7201          FORMAT(T3,' MANR scale factor for NR step is ',
     &                F5.3,' for ', A5,'.')
               SCRATCH(NOPT + IOPT) = SCRATCH(NOPT + IOPT)*CORR
            ENDIF
C     
         ENDIF
C
      ENDDO 
C
      RETURN
      END


