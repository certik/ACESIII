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

C THIS PROJECTS OUT LINEARLY DEPENDENT PARTS OF THE INTERNAL COORDINATES

c INPUT
c integer NATOM

c OUTPUT
c double SCR(*)
c double SYMQ(*)

c RECORDS
c get 'COMPNSYQ'
c get 'COMPSYMQ'
c put 'COMPSYMQ'

      SUBROUTINE SCHMIDT(NATOM,SCR,SYMQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SCR(*),SYMQ(*)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA TOL /1.D-8/

      NSIZE=3*NATOM

      CALL IGETREC(20,'JOBARC','COMPNSYQ',1,NSYMOLD)
      CALL DGETREC(20,'JOBARC','COMPSYMQ',NSIZE*NSIZE,SCR)

c   o loop over vectors
      DO I=2,NSYMOLD
         IPOS=(I-1)*NSIZE+1
         CALL GSCHMIDT(SCR(IPOS),SCR,NSIZE,I-1,SYMQ,ZJUNK,1.D-8)
      END DO
      CALL XDCOPY(NSIZE*NSIZE,SCR,1,SYMQ,1)

c   o renormalize
      IOFF=1
      DO I=1,NSYMOLD
         X=DNRM2(NSIZE,SYMQ(IOFF),1)
         IF (ABS(X).GT.TOL) THEN
            Z=1.d0/X
            CALL XDSCAL(NSIZE,Z,SYMQ(IOFF),1)
         END IF
         IOFF=IOFF+NSIZE
      END DO

      CALL DPUTREC(20,'JOBARC','COMPSYMQ',NSIZE*NSIZE,SYMQ)

      RETURN
      END

