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
      SUBROUTINE PERPOP(T,QREF,SCRATCH,ATMASS,NATOM,NORD,SYMTOL,IFOUND)
C
C THIS ROUTINE CHECKS TO SEE IF THERE IS A C2 AXIS IN THE XY
C PLANE FOR THE COORDINATES GIVEN IN QREF.  
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 T
      DIMENSION QREF(3*NATOM),ATMASS(NATOM),NORD(NATOM),SCRATCH(1)
      DIMENSION RM(9)
C
      DATA ONE,ZILCH /1.0D+00,0.0D+00/
C
      IFOUND=0
      I1=3*NATOM+1
C
C FIRST IDENTIFY AN ATOM WHICH HAS A NONZERO PROJECTION IN THE XY PLANE
C
      IPERM=0
      DO 10 I=1,NATOM
       XCOR=QREF((I-1)*3+1)*ATMASS(I)
       YCOR=QREF((I-1)*3+2)*ATMASS(I)
       ZCOR=QREF((I-1)*3+3)*ATMASS(I)
       IF(MAX(ABS(XCOR),ABS(YCOR)).GT.SYMTOL.AND.IPERM.EQ.0)THEN
        IPERM=I
        ZMATCH=ZCOR
        IPOS0=(I-1)*3+1
       ENDIF
10    CONTINUE
C
C NOW LOOP OVER ALL ATOMS HAVING THE SAME ATOMIC NUMBER AND Z
C  COORDINATES WITH OPPOSITE SIGNS AND TEST FOR C2 THROUGH BISECTOR
C
      DO 11 I=1,NATOM
       IF(I.NE.IPERM.AND.ATMASS(I).NE.ZILCH)THEN
        ZCOR=QREF((I-1)*3+3)*ATMASS(I)
        IF(T.EQ.'C')CHECK=ZCOR+ZMATCH
        IF(T.EQ.'P')CHECK=ZCOR-ZMATCH
        IF(ABS(CHECK).LT.SYMTOL)THEN
C
C FIND POSITION OF BISECTING LINE AND ROTATE SYSTEM SO THAT IT LIES
C  ALONG THE X AXIS.  
C
         CALL ALONGX(QREF((I-1)*3+1),QREF(IPOS0),RM,IERR)
         IF(IERR.NE.-1)THEN
          CALL MTRAN2(RM,3)
C
CJDW 1/6/98. Replace MATMULV call by XGEMM call. SCRATCH=(RM)^T * QREF.
C
C         CALL MATMULV(SCRATCH,QREF,RM,NATOM,3,3)
C
          CALL XGEMM('T','N',3,NATOM,3,ONE,RM,3,QREF,3,ZILCH,SCRATCH,3)
C
C NOW CHECK TO SEE IF THIS IS A AXIS OR AN X->-X PLANE
C
          IF(T.EQ.'C')THEN
           CALL DOSYOP('C',2,1,1,NATOM,SCRATCH,SCRATCH(I1),IDUM,0)
          ELSEIF(T.EQ.'P')THEN
           CALL DOSYOP('P',1,1,2,NATOM,SCRATCH,SCRATCH(I1),IDUM,0)
          ELSE
           WRITE(6,1000)T
1000       FORMAT(T3,'@PERPOP-F, Operation type ',A,' unknown.')
           CALL ERREX
          ENDIF
          CALL COMPARE(SCRATCH(I1),SCRATCH,ATMASS,NORD,NATOM,ICOMPX,
     &                  SYMTOL)
          IF(ICOMPX.EQ.0)THEN
           CALL XDCOPY(3*NATOM,SCRATCH,1,QREF,1)
           IFOUND=1
          ENDIF
         ENDIF
        ENDIF
       ENDIF
11    CONTINUE
C
      RETURN
      END   
