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

      SUBROUTINE POLYPRT0(NATOMS,IINTFP,IPRINT,
     &                    LCOORD,LENERG,LGRAD,LHESS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LCOORD,LENERG,LGRAD,LHESS
      PARAMETER(MXATMS=250)
      DIMENSION  SCR(3*MXATMS + 3*MXATMS + 9*MXATMS*MXATMS)
      DIMENSION ISCR(  MXATMS)

      IF (NATOMS.GT.MXATMS) THEN
         WRITE(6,1000) NATOMS, MXATMS
 1000    FORMAT('  @POLYPRT0: Too many atoms. Given ',I3,' Max ',I3)
         CALL ERREX
      END IF

      I000 = 1
      I010 = I000 + 3*NATOMS
      I020 = I010 + 3*NATOMS
      I030 = I020 + 9*NATOMS*NATOMS
      CALL POLYPRT(SCR(I000),SCR(I010),SCR(I020),ISCR,
     &             NATOMS,IINTFP,IPRINT,
     &             LCOORD,LENERG,LGRAD,LHESS)

      RETURN
      END
