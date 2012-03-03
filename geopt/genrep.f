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





      SUBROUTINE GENREP(BIGX,R,IPTR,NATOM)
C
C GENERATES THE (3*NATOM,3*NATOM) REPRESENTATION OF
C THE SYMMETRY OPERATION GIVEN BY THE 3x3 CARTESIAN
C MATRIX R AND THE POINTER VECTOR IPTR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(3,3),IPTR(NATOM),BIGX(3*NATOM,3*NATOM)
      DATA ONE,ONEM,ZILCH/1.0D0,-1.0D0,0.0D0/
      CALL ZERO(BIGX,9*NATOM*NATOM)
      NX=3*NATOM
c      CALL MTRAN2(R,3)
      DO 10 IATOM=1,NATOM
       JATOM=IPTR(IATOM)
       I=1+3*(IATOM-1)
       J=1+3*(JATOM-1)
       CALL BLKCPY(R,3,3,BIGX,NX,NX,I,J)
10    CONTINUE
      RETURN
      END  
