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
      SUBROUTINE COMPARE(Q1,Q2,ATMASS,IMAP,NATOM,IERR,TOL)
C
C THIS ROUTINE CHECKS WHETHER COORDINATE VECTORS Q1 AND
C Q2 ARE RELATED BY A SIMPLE PERMUTATION.  IF THEY ARE,
C IERR IS RETURNED AS "0", AND IMAP CONTAINS THE PERMUTATION
C VECTOR
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IMATCH
      DIMENSION Q1(3*NATOM),Q2(3*NATOM),ATMASS(NATOM),IMAP(NATOM)
      DIMENSION TMP(3)
C
      DATA ZILCH / 0.0/
      DATA ONEM  /-1.0/ 
C
      IERR=0
      IOFF1=1
c      CALL IZERO(IMAP,NATOM)
      do i = 1, natom
         imap(i) = 0
      enddo

      DO 10 IATOM1=1,NATOM
       IF(ATMASS(IATOM1).NE.ZILCH)THEN
        IOFF2=1
        IMATCH=.FALSE.
        DO 20 IATOM2=1,NATOM
         IF(ATMASS(IATOM2).EQ.ATMASS(IATOM1))THEN
          CALL VADD(TMP,Q1(IOFF1),Q2(IOFF2),3,ONEM)
          I=ISAMAX(3,TMP,1)
          IF(ABS(TMP(I)).LT.TOL)THEN
           IMAP(IATOM1)=IATOM2
           IMATCH=.TRUE.
          ENDIF
         ENDIF
         IOFF2=IOFF2+3
20      CONTINUE
        IF(.NOT.IMATCH)IERR=IERR+1
       ENDIF
       IOFF1=IOFF1+3
10    CONTINUE
C
      RETURN
      END
