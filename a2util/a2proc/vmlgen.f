C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      SUBROUTINE VMLGEN(PGRP,IORDGP,QREF,QGEN,IGEN)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: vmlgen.f,v 1.2 2006/02/24 00:16:21 yau Exp $
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
       CALL SCOPY(3,QREF,1,QGEN(1,I),1)
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
       CALL GETREC(20,'JOBARC','D2TYPXYZ',1,ITYPE)
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
       CALL SCOPY(3,QGEN(1,I),1,QGEN(1,IGEN),1)
20    CONTINUE
      RETURN
      END
