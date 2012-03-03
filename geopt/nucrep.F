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
      SUBROUTINE NUCREP(ZNREP,IERR)
C
C DETERMINES THE INTERNUCLEAR REPULSION ENERGY.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Maximum number of atoms currently allowed
#include "mxatms.par"
#include "coord.com"
C
      IQ(I)=3*I-2
      Z=0.D0
      DO 10 I=1,NATOMS
      DO 10 J=I+1,NATOMS
       BOT=DIST(Q(IQ(I)),Q(IQ(J)))
       IF(BOT.LT.1.D-5)THEN
        IERR=1
        GOTO 10
       ELSE
        IF(IATNUM(I).NE.110.AND.IATNUM(J).NE.110)THEN
c        IGNORE GHOST ATOMS
         TOP=FLOAT(IATNUM(I)*IATNUM(J))
         Z=Z+(TOP/BOT)
        ENDIF
      ENDIF
10    CONTINUE
      ZNREP=Z
      RETURN
      END
