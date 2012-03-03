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
      SUBROUTINE FRMLST(NATOM,IREF,IORBIT,IORDGP,TRNMAT,IATLST,IPTR,
     &                  MEMBER,ICOUNT,JCOUNT)
C
C THIS ROUTINE CREATES THE MEMBER AND TRNMAT ARRAYS FOR A PARTICULAR
C  ORBIT, GIVEN THE REFERENCE ATOM AND IPTR ARRAYS. 
C
      IMPLICIT INTEGER (A-H,O-Z)
      DIMENSION IPTR(NATOM,IORDGP),IATLST(NATOM) 
      INTEGER TRNMAT(NATOM),MEMBER(NATOM)
      PARAMETER (LUOUT = 6)
      JIN=0  
      IFIRST=0
      DO 10 I=1,IORDGP
       IATLST(IPTR(IREF,I))=IORBIT
       TRNMAT(IPTR(IREF,I))=I
       IF(IREF.EQ.IPTR(IREF,I))TRNMAT(IPTR(IREF,I))=0
10    CONTINUE
      DO 33 I=1,NATOM
       IF(IATLST(I).EQ.IORBIT)THEN
        JIN=JIN+1
        ICOUNT=ICOUNT+1
        JCOUNT=JCOUNT+1
        IF(JIN.EQ.1)IFIRST=ICOUNT
        MEMBER(ICOUNT)=I
       ENDIF
 33   CONTINUE
C
C MAKE SURE THAT THE FIRST ENTRY FOR THIS ORBIT IN MEMBER IS THE
C  REFERENCE ATOM. IF NOT, SWAP THE ORDERING.
C
      JLOC=CHKINT(MEMBER,NATOM,IREF)
      IF(JLOC.NE.IFIRST)THEN
       WRITE(LUOUT,100)
100    FORMAT(T3,'@FRMLST-I, Member ordering differs from Z-matrix',
     &           ' order.')
       ITEMP=MEMBER(IFIRST)
       MEMBER(IFIRST)=IREF
       MEMBER(JLOC)=ITEMP
      ENDIF
      RETURN
      END
