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
      SUBROUTINE VMLGENX(PGRP,IORDGP,QREF,QGEN,IGEN, ITYPE)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: vmlgen.f,v 1.3 2010/06/30 15:58:42 ponton Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     vmlgen -- generate symmetry equivalent atoms like VMol does
C
C SYNOPSIS
      CHARACTER*4 PGRP
      Integer IOrdGp, IGen, ITYPE
      Double precision QREF(3),QGEN(3,IORDGP-1)
C
C DESCRIPTION
C     Given one set of cartesian coordinates produces all symmetry
C     equivalent coordinates in the same fashion as VMol uses internally
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer I, j
      LOGICAL QSAME
      DO 10 I=1,IORDGP-1
        QGEN(1,I) = QREF(1) 
        QGEN(2,I) = QREF(2) 
        QGEN(3,I) = QREF(3) 
10    CONTINUE
      IF(PGRP(1:3).EQ.'C2v')THEN
       QGEN(1,1)=-QGEN(1,1)
       QGEN(2,2)=-QGEN(2,2)
       QGEN(1,3)=-QGEN(1,3)
       QGEN(2,3)=-QGEN(2,3)
      ELSEIF(PGRP(1:3).EQ.'C2h')THEN
       QGEN(3,1)=-QGEN(3,1)
       QGEN(1,2)=-QGEN(1,2)
       QGEN(2,2)=-QGEN(2,2)
       QGEN(1,3)=-QGEN(1,3)
       QGEN(2,3)=-QGEN(2,3)
       QGEN(3,3)=-QGEN(3,3)
      ELSEIF(PGRP(1:3).EQ.'D2 ')THEN 
C
C What was here was inadequate and I believe that it has been like
C that for a long time. Thankfully, we do not have a lot of D2
C molecules. 01/2006, Ajith Perera.
C
       IF (ITYPE.EQ.1) THEN
          QGEN(1,1)=-QGEN(1,1)
          QGEN(2,1)=-QGEN(2,1)
          QGEN(1,2)=-QGEN(1,2)
          QGEN(3,2)=-QGEN(3,2)
          QGEN(2,3)=-QGEN(2,3)
          QGEN(3,3)=-QGEN(3,3)
       ELSE IF (ITYPE.EQ.2) THEN
          QGEN(1,1)=-QGEN(1,1)
          QGEN(2,1)=-QGEN(2,1)
          QGEN(2,2)=-QGEN(2,2)
          QGEN(3,2)=-QGEN(3,2)
          QGEN(1,3)=-QGEN(1,3)
          QGEN(3,3)=-QGEN(3,3)
       ELSE IF (ITYPE.EQ.3) THEN
          QGEN(1,1)=-QGEN(1,1)
          QGEN(3,1)=-QGEN(3,1)
          QGEN(2,2)=-QGEN(2,2)
          QGEN(3,2)=-QGEN(3,2)
          QGEN(1,3)=-QGEN(1,3)
          QGEN(2,3)=-QGEN(2,3)
       END IF
C
      ELSEIF(PGRP(1:3).EQ.'C2 ')THEN
       QGEN(1,1)=-QGEN(1,1)
       QGEN(2,1)=-QGEN(2,1)
      ELSEIF(PGRP(1:3).EQ.'C s')THEN
       QGEN(3,1)=-QGEN(3,1)
      ELSEIF(PGRP(1:3).EQ.'C i')THEN
       QGEN(1,1)=-QGEN(1,1)
       QGEN(2,1)=-QGEN(2,1)
       QGEN(3,1)=-QGEN(3,1)
      ELSEIF(PGRP(1:3).EQ.'D2h')THEN
C
C I checked this. It is correct. 01/2006, Ajith Perera.
C
       QGEN(1,1)=-QGEN(1,1)
       QGEN(2,2)=-QGEN(2,2)
       QGEN(1,3)=-QGEN(1,3)
       QGEN(2,3)=-QGEN(2,3)
       QGEN(3,4)=-QGEN(3,4)
       QGEN(1,5)=-QGEN(1,5)
       QGEN(3,5)=-QGEN(3,5)
       QGEN(2,6)=-QGEN(2,6)
       QGEN(3,6)=-QGEN(3,6)
       QGEN(1,7)=-QGEN(1,7)
       QGEN(2,7)=-QGEN(2,7)
       QGEN(3,7)=-QGEN(3,7)
      ENDIF
C
C NOW NOT ALL OF THESE CENTERS WILL BE UNIQUE.  DETERMINE WHICH ONES
C  ARE UNIQUE AND STUFF THESE POSITIONS AT THE BOTTOM OF QGEN
C
      IGEN=0 
      DO 20 I=1,IORDGP-1
C
C SEE IF THE GENERATED CENTER IS EQUAL TO THE REFERENCE CENTER
C
       IF(QSAME(QREF,QGEN(1,I)))GOTO 20
       DO 22 J=1,I-1
C
C SEE IF THE GENERATED CENTER IS EQUAL TO A PREVIOUSLY GENERATED
C   REDUNDANT CENTER
C
        IF(QSAME(QGEN(1,J),QGEN(1,I)))GOTO 20
22     CONTINUE
C
C THIS IS A UNIQUE CENTER.  COPY IT TO THE FIRST AVAILABLE POSITION IN QGEN
C
       IGEN=IGEN+1
       QGEN(1,IGEN)=QGEN(1,I)
       QGEN(2,IGEN)=QGEN(2,I)
       QGEN(3,IGEN)=QGEN(3,I)
20    CONTINUE
      RETURN
      END
