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
       SUBROUTINE GENRST(J,SCR,NTIME,NUMTIM,IORDGP,JX,CLSTYP)
C
C THIS SUBROUTINE GENERATES THE REMAINING TRANSFORMATION MATRICES
C  BELONGING TO A GIVEN CLASS BY SUCCESSIVELY APPLYING A PARTICULAR
C  TRANSFORMATION MATRIX.  PARTICULARLY USEFUL FOR NON-CUBIC GROUPS.
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION SCR(9*IORDGP)
       INTEGER CLSTYP(IORDGP) 
       IP(I)=1+9*(I-1) 
       DO 10 I=1,NUMTIM
        NTIME=NTIME+1
        CALL UNITRY(SCR(IP(J)),SCR(IP(NTIME-1)),SCR(IP(NTIME)),3)
        CLSTYP(NTIME)=JX
10     CONTINUE
       RETURN
       END
