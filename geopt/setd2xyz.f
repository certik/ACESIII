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
      SUBROUTINE SETD2XYZ(NATOM,Q,ATMASS,STRING1,STRING2,IERROR)
C
C SELECT PROPER C2 AXES TO WRITE TO MOL FILE.  THIS ROUTINE IS
C NECESSARY BECAUSE VMOL CRASHES WHEN THE TWO OPERATIONS IN THE
C MOL FILE MAP ANY OF THE ATOMS TO THE SAME LOCATION.  THIS IS
C NOT A PROBLEM FOR C2H AND C2V, BUT CAN CAUSE SEVERE FRUSTRATION
C FOR D2.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Q(3,NATOM),ATMASS(NATOM)
      LOGICAL ATOMONX,ATOMONY,ATOMONZ,XATOM,YATOM,ZATOM
      CHARACTER*3 STRING1,STRING2
      PARAMETER (TOL=1.D-6)
      ATOMONX=.FALSE.
      ATOMONY=.FALSE.
      ATOMONZ=.FALSE.
      IERROR=0
      DO 10 IATOM=1,NATOM
       IF(ATMASS(IATOM).NE.0.0D0)THEN
        XATOM=DABS(Q(1,IATOM)).GT.TOL
        YATOM=DABS(Q(2,IATOM)).GT.TOL
        ZATOM=DABS(Q(3,IATOM)).GT.TOL
        IF(XATOM.AND..NOT.YATOM.AND..NOT.ZATOM)ATOMONX=.TRUE.
        IF(YATOM.AND..NOT.XATOM.AND..NOT.ZATOM)ATOMONY=.TRUE.
        IF(ZATOM.AND..NOT.YATOM.AND..NOT.XATOM)ATOMONZ=.TRUE.
       ENDIF
10    CONTINUE
      IF(.NOT.ATOMONX)THEN
       STRING1=' XY'
       STRING2=' XZ'
       ITYPE = 1
      ELSEIF(.NOT.ATOMONY)THEN
       STRING1=' XY'
       STRING2=' YZ'
       ITYPE = 2
      ELSEIF(.NOT.ATOMONZ)THEN
       STRING1=' XZ'
       STRING2=' YZ'
       ITYPE = 3
      ELSE
       WRITE(6,1000)
       IERROR=1   
       ITYPE = 0
       CALL ERREX
      ENDIF
C
C We need to know the type of D2 axis to make the correct map from vmol
C to zmat in vmol2ja. Otherwise all of the geometry optimizations and
C Freq. for D2 point group are going to suffer. 01/2006, Ajith Perera.
C
      CALL IPUTREC(20, 'JOBARC', 'D2TYPXYZ', 1, ITYPE)
C
      RETURN
1000  FORMAT(T3,'@SETD2XYZ-I, The integral program is unable to use ',
     &          'D2 symmetry for this case.')
      END
