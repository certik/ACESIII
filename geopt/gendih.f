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
      SUBROUTINE GENDIH(OPRTYP,OPRORD,OPRTIM,OPRAXS,NUMCLS,
     &                  NORDER,NTIME,TRNMAT)
C
C THIS ROUTINE GENERATES TRANSFORMATION MATRICES FOR DIHEDRAL
C  PLANES AND DIHEDRAL C2 ROTATION AXES.  THESE ARE NOT COINCIDENT
C  WITH THE CARTESIAN AXIS SYSTEM OF THE CANONICAL JODA ORIENTATION
C  AND HENCE ARE A PAIN IN THE ASS.  
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TRNMAT(3,3),SCR1(3,3),SCR2(3,3)
      INTEGER OPRORD(2),OPRTIM(2),OPRAXS(2),NUMCLS(2)
      CHARACTER*1 OPRTYP(2)
      IERR=0
C
C GET TRANS. MAT. FOR THE V, RATHER THAN THE H OPERATION.
C
      CALL IDNMAT(SCR1,3,NTIME)
      CALL DOSYOP(OPRTYP(1),OPRORD(1),OPRTIM(1),OPRAXS(1),
     &            3,SCR1,TRNMAT,IERR,1)
C
C GET APPROPRIATE ROTATION MATRIX FOR TRANSFORMATION.
C
      NTIME=NTIME-1
      CALL IDNMAT(SCR1,3,NTIME)
      CALL DOSYOP('C',2*NORDER,1,3,3,SCR1,SCR2,IERR,1)
C
C TRANSFORM THE V OP TO A D OP.
C
      CALL UNITRY(SCR2,TRNMAT,TRNMAT,3)
C
C NOW FILL IN THE CLASS VECTORS WITH INFO ABOUT THIS OP.
C
      OPRTYP(2)=OPRTYP(1)
      OPRORD(2)=OPRORD(1)
      OPRAXS(2)=4
      OPRTIM(2)=OPRTIM(1)
      NUMCLS(2)=NUMCLS(1)
      RETURN
      END
