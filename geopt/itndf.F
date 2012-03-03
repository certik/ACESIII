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
C
C DETERMINES THE TOTAL NUMBER OF DEGREES OF FREEDOM WITHIN
C  THE TOTALLY SYMMETRIC MOLECULAR SUBSPACE.  EQUAL TO THE
C  NUMBER OF TOTALLY SYMMETRIC VIBRATIONS.
C
C ALGORITHM LOOSELY BASED ON A PAPER BY FOWLER AND QUINN
C   IN THEO. CHEM. ACTA (1986).
C
      INTEGER FUNCTION ITNDF()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
      CHARACTER*8 STSYM(MXATMS)
      character*(8*mxatms) szStSymTmp
      LOGICAL YESNO
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C PULL NUMBER OF ORBITS AND SITE GROUP VECTOR FROM JOBARC.
C
      ISIZE=1
      CALL IGETREC(20,'JOBARC','FULLNORB',ISIZE,IORBIT)
      CALL GETCREC(20,'JOBARC','FULLSTGP',IORBIT*8,szStSymTmp)
c   o 'compress' szStSymTmp into the 2-D STSYM array
      iNdx = 1
      do i = 1, IORBIT
         STSYM(i)(1:8) = szStSymTmp(iNdx:iNdx+7)
         iNdx = iNdx + 8
      end do
C
C COMPUTE THE NUMBER OF TRANSLATIONS AND ROTATIONS WHICH
C  TRANSFORM AS THE TOTALLY SYMMETRIC IR.
C
      IVR=NRTSRV(FPGRP)
C
C GO THROUGH IT AND ASSIGN THE PROPER NUMBER OF D.O.F.
C  PER ORBIT.
C
      ITNDF=0
      YESNO=.FALSE.
      DO 10 I=1,IORBIT
       INCREM=1
       IF(STSYM(I)(1:3).EQ.'C1 ')INCREM=3
       IF(STSYM(I)(1:3).EQ.'C s')INCREM=2
       IF((FPGRP(1:1).NE.'C'.OR.
     &    FPGRP.EQ.'C i').AND.STSYM(I)(1:4).EQ.FPGRP)INCREM=0
       ITNDF=ITNDF+INCREM
10    CONTINUE
      ITNDF=ITNDF-IVR

      RETURN
      END
