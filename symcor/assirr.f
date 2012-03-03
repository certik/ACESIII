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
      SUBROUTINE ASSIRR(PTGRPX,CHRVEC,NDEG,IRREP)
C
C THIS FUNCTION RETURNS THE IRREDUCIBLE REPRESENTATION CORRESPONDING
C TO A VECTOR OF LENGTH H (THE ORDER OF THE GROUP) WHICH CONTAINS THE
C CHARACTERS UNDER ALL OPERATIONS.
C
C INPUT:
C    PTGRPX (CHARACTER*4) - THE POINT GROUP OF THE MOLECULE
C    CHRVEC               - A VECTOR WHICH CONTAINS THE CHARACTERS
C                           OF ALL OPERATIONS IN THE POINT GROUP
C
C RETURNED:
C
C    IRREP                - THE IRREDUCIBLE REPRESENTATION OF V
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ATOI
      DIMENSION ANGLE(50)
      PARAMETER (MXATMS=50)
      PARAMETER (TOL = 1.D-2)
      CHARACTER*4 PTGRP,PTGRPX,IRREP
      CHARACTER*2 GREEK(6)     
      DIMENSION CHRVEC(*)
      DATA GREEK /'SG','PI','DE','PH','GA','IO'/
C
C SET INTEGER TO POSITION OF LAST NONBLANK CHARACTER IN THE STRING
C
      ICLST=linblnk(PTGRPX)
C
C SOME USEFUL CONSTANTS
C
      PI=DACOS(-1.D0)
      TWOPI=2.D0*PI
C
C INITIALIZE IRREP TO BLANKS
C
      IRREP='    '
C
C**********************************************************************
C DEAL FIRST WITH NON-CUBIC GROUPS WHICH CAN NOT BE HANDLED BY
C THE GENERAL C(nh), C(nv), D(nh), D(nd), D(n), C(n) AND S(2n)
C ALGORITHMS.  TREAT Cs AS C1h, Ci AS S2.
C**********************************************************************
C
      PTGRP=PTGRPX
      IF(PTGRPX(1:3).EQ.'C s')PTGRP='C1h '
      IF(PTGRPX(1:3).EQ.'C i')PTGRP='S2  '
C
C**********************************************************************
C WE REQUIRE SPECIAL CODE FOR D2 AND D2h.  HERE IT IS.
C**********************************************************************
C
      IF(PTGRP(1:3).EQ.'D2 ')THEN
       CHARC2Z=CHRVEC(1)
       CHARC2X=CHRVEC(2)
       IF(DABS(CHARC2Z-1.D0).LT.TOL.AND.DABS(CHARC2X-1.D0).LT.TOL)THEN
        IRREP(4:4)='A'
       ENDIF
       IF(DABS(CHARC2Z-1.D0).LT.TOL.AND.DABS(CHARC2X+1.D0).LT.TOL)THEN
        IRREP(3:4)='B1'
       ENDIF
       IF(DABS(CHARC2Z+1.D0).LT.TOL.AND.DABS(CHARC2X+1.D0).LT.TOL)THEN
        IRREP(3:4)='B2'
       ENDIF
       IF(DABS(CHARC2Z+1.D0).LT.TOL.AND.DABS(CHARC2X-1.D0).LT.TOL)THEN
        IRREP(3:4)='B3'
       ENDIF
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'D2h')THEN
       CHARC2Z=CHRVEC(1)
CJDW 8/4/97. I think next line is wrong.
C      CHARC2X=CHRVEC(4)
C            I think it should be :
       CHARC2X=CHRVEC(3)
       CHARINV=CHRVEC(2)
       IF(DABS(CHARINV-1.D0).LT.TOL)IRREP(4:4)='g'
       IF(DABS(CHARINV+1.D0).LT.TOL)IRREP(4:4)='u'
       IF(DABS(CHARC2Z-1.D0).LT.TOL.AND.DABS(CHARC2X-1.D0).LT.TOL)THEN
        IRREP(3:3)='A'
       ENDIF
       IF(DABS(CHARC2Z-1.D0).LT.TOL.AND.DABS(CHARC2X+1.D0).LT.TOL)THEN
        IRREP(2:3)='B1'
       ENDIF
       IF(DABS(CHARC2Z+1.D0).LT.TOL.AND.DABS(CHARC2X+1.D0).LT.TOL)THEN
        IRREP(2:3)='B2'
       ENDIF
       IF(DABS(CHARC2Z+1.D0).LT.TOL.AND.DABS(CHARC2X-1.D0).LT.TOL)THEN
        IRREP(2:3)='B3'
       ENDIF
       RETURN
C
C**********************************************************************
C  THE FOLLOWING CODE PERTAINS TO THE INFINITE GROUPS CXv AND DXh.
C**********************************************************************
C
      ELSEIF(PTGRP(1:3).EQ.'DXh')THEN
       CHAR=CHRVEC(14)
       CHAR2=CHRVEC(31)
       CHAR3=CHRVEC(1)
       IF(NDEG.EQ.1)IRREP(1:2)='SG'
       IF(NDEG.EQ.2)THEN
       IRREP(1:2)='PI'
        DO 3100 J=1,4
         Z=2.D0*DCOS(DFLOAT(J)*DACOS(-1.D0)/4.D0)
         IF(DABS(CHAR3-Z).LT.TOL)IRREP(1:2)=GREEK(J+1)
 3100   CONTINUE
        CHAR=CHAR/2.0
       ENDIF
       IF(DABS(CHAR+1.D0).LT.TOL)IRREP(3:3)='u'
       IF(DABS(CHAR-1.D0).LT.TOL)IRREP(3:3)='g'
       IF(NDEG.EQ.1)THEN
        IF(DABS(CHAR+1.D0).LT.TOL)IRREP(4:4)='-'
        IF(DABS(CHAR-1.D0).LT.TOL)IRREP(4:4)='+'
       ENDIF
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'CXv')THEN
       CHAR=CHRVEC(15)
       CHAR3=CHRVEC(1)
       IF(NDEG.EQ.1)THEN
        IRREP(2:3)='SG'
        IF(DABS(CHAR+1.D0).LT.TOL)IRREP(4:4)='-'
        IF(DABS(CHAR-1.D0).LT.TOL)IRREP(4:4)='+'
       ELSE
        IRREP(3:4)='PI'
        DO 3101 J=1,3
         Z=2.D0*DCOS(DFLOAT(J)*DACOS(-1.D0)/4.D0)
         IF(DABS(CHAR3-Z).LT.TOL)IRREP(3:4)=GREEK(J+1)
 3101   CONTINUE
       ENDIF
C
C**********************************************************************
C C(nh) GROUPS.  SEPARATE CODE FOR n ODD AND EVEN
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'C'.AND.PTGRP(ICLST:ICLST).EQ.'h')THEN
       IORDER=ATOI(PTGRP(2:ICLST-1))
       ITOP=(IORDER-1)/2
       DO 100 I=1,ITOP
        ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
100    CONTINUE
       IF(MOD(IORDER,2).NE.0)THEN
        CHARCN=CHRVEC(1)
        CHARSGH=CHRVEC(IORDER*2)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARSGH=CHARSGH/DFLOAT(NDEG)
        IF(DABS(CHARSGH-1.0D0).LT.TOL)IRREP(3:4)=''''
        IF(DABS(CHARSGH+1.0D0).LT.TOL)IRREP(3:4)=''''''
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARCN-1.D0).LT.TOL)IRREP(2:2)='A'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(1:1)='E'
          DO 101 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(2:2),'(I1)')I
101       CONTINUE
         ELSE
          IRREP(2:2)='E'
         ENDIF
        ENDIF
C
       ELSE 
C
        CHARCN=CHRVEC(1)
        CHARINV=CHRVEC(IORDER+IORDER/2)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARINV=CHARINV/DFLOAT(NDEG)
        IF(DABS(CHARINV-1.D0).LT.TOL)IRREP(4:4)='g'
        IF(DABS(CHARINV+1.D0).LT.TOL)IRREP(4:4)='u'
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARCN-1.D0).LT.TOL)IRREP(3:3)='A'
         IF(DABS(CHARCN+1.D0).LT.TOL)IRREP(3:3)='B'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(2:2)='E'
          DO 102 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(3:3),'(I1)')I
102       CONTINUE
         ELSE
          IRREP(3:3)='E'
         ENDIF
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C C(nv) GROUPS.  SEPARATE CODE FOR n ODD AND EVEN
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'C'.AND.PTGRP(ICLST:ICLST).EQ.'v')THEN
       IORDER=ATOI(PTGRP(2:ICLST-1))
       ITOP=(IORDER-1)/2
       DO 200 I=1,ITOP
        ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
200    CONTINUE
       IF(MOD(IORDER,2).NE.0)THEN
        CHARCN=CHRVEC(1)
        CHARSGV=CHRVEC(IORDER)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARSGV=CHARSGV/DFLOAT(NDEG)
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARSGV-1.D0).LT.TOL)IRREP(3:4)='A1'
         IF(DABS(CHARSGV+1.D0).LT.TOL)IRREP(3:4)='A2'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(3:3)='E'
          DO 201 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(4:4),'(I1)')I
201       CONTINUE
         ELSE
          IRREP(4:4)='E'
         ENDIF
        ENDIF
C
       ELSE 
C
        CHARCN=CHRVEC(1)
        CHARSGD=CHRVEC(IORDER+(IORDER/2))
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARSGD=CHARSGD/DFLOAT(NDEG)
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARCN-1.D0).LT.TOL)IRREP(3:3)='A'
         IF(DABS(CHARCN+1.D0).LT.TOL)IRREP(3:3)='B'
         IF(DABS(CHARSGD-1.D0).LT.TOL)IRREP(4:4)='1'
         IF(DABS(CHARSGD+1.D0).LT.TOL)IRREP(4:4)='2'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(3:3)='E'
          DO 202 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(4:4),'(I1)')I
202       CONTINUE
         ELSE
          IRREP(4:4)='E'
         ENDIF
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C D(nd) GROUPS.  SEPARATE CODE FOR n ODD AND EVEN
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'D'.AND.PTGRP(ICLST:ICLST).EQ.'d')THEN
       IORDER=ATOI(PTGRP(2:ICLST-1))
       IF(MOD(IORDER,2).NE.0)THEN
        ITOP=(IORDER-1)/2
        DO 300 I=1,ITOP
         ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
300     CONTINUE
        CHARCN=CHRVEC(1)
        CHARC2=CHRVEC(2*IORDER)
        CHARINV=CHRVEC(2*IORDER-1)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARC2=CHARC2/DFLOAT(NDEG)
        CHARINV=CHARINV/DFLOAT(NDEG)
        IF(DABS(CHARINV-1.0D0).LT.TOL)IRREP(4:4)='g'
        IF(DABS(CHARINV+1.0D0).LT.TOL)IRREP(4:4)='u'
        IF(NDEG.EQ.1)THEN
         IRREP(2:2)='A'
         IF(DABS(CHARC2-1.D0).LT.TOL)IRREP(3:3)='1'
         IF(DABS(CHARC2+1.D0).LT.TOL)IRREP(3:3)='2'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(2:2)='E'
          DO 301 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(3:3),'(I1)')I
301       CONTINUE
         ELSE
          IRREP(3:3)='E'
         ENDIF
        ENDIF
C
       ELSE 
C
        ITOP=(2*IORDER-1)/2
        DO 302 I=1,ITOP
         ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(2*IORDER)
302     CONTINUE
        CHARCN=CHRVEC(1)
        CHARS2N=CHRVEC(IORDER)
        CHARC2=CHRVEC(2*IORDER)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARC2=CHARC2/DFLOAT(NDEG)
        CHARS2N=CHARS2N/DFLOAT(NDEG)
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARS2N-1.D0).LT.TOL)IRREP(3:3)='A'
         IF(DABS(CHARS2N+1.D0).LT.TOL)IRREP(3:3)='B'
         IF(DABS(CHARC2-1.D0).LT.TOL)IRREP(4:4)='1'
         IF(DABS(CHARC2+1.D0).LT.TOL)IRREP(4:4)='2'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN 
          DO 303 I=1,ITOP
           IRREP(3:3)='E'
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARS2N-Z).LT.TOL)WRITE(IRREP(4:4),'(I1)')I
303       CONTINUE
         ELSE
          IRREP(4:4)='E'
         ENDIF
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C D(nh) GROUPS.  SEPARATE CODE FOR n ODD AND EVEN. 
C   DOES NOT WORK FOR D2h
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'D'.AND.PTGRP(ICLST:ICLST).EQ.'h')THEN
       IORDER=ATOI(PTGRP(2:ICLST-1))
       ITOP=(IORDER-1)/2
       DO 400 I=1,ITOP
        ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
400     CONTINUE
       IF(MOD(IORDER,2).NE.0)THEN
        CHARCN=CHRVEC(1)
        CHARC2=CHRVEC(2*IORDER-1)
        CHARSGH=CHRVEC(4*IORDER-1)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARC2=CHARC2/DFLOAT(NDEG)
        CHARSGH=CHARSGH/DFLOAT(NDEG)
        IF(NDEG.EQ.1)THEN
         IRREP(1:1)='A'
         IF(DABS(CHARC2-1.D0).LT.TOL)IRREP(2:2)='1'
         IF(DABS(CHARC2+1.D0).LT.TOL)IRREP(2:2)='2'
         IF(DABS(CHARSGH-1.0D0).LT.TOL)IRREP(3:4)=''''
         IF(DABS(CHARSGH+1.0D0).LT.TOL)IRREP(3:4)=''''''
        ELSEIF(NDEG.EQ.2)THEN
         IF(DABS(CHARSGH-1.0D0).LT.TOL)IRREP(3:4)=''''
         IF(DABS(CHARSGH+1.0D0).LT.TOL)IRREP(3:4)=''''''
         IF(ITOP.NE.1)THEN
          IRREP(1:1)='E'
          DO 401 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(2:2),'(I1)')I
401       CONTINUE
         ELSE
          IRREP(2:2)='E'
         ENDIF
        ENDIF
C
       ELSE 
C
        CHARCN=CHRVEC(1)
        CHARC2=CHRVEC(2*IORDER-1)
        CHARINV=CHRVEC(2*IORDER-2)
        CHARCN=CHARCN/DFLOAT(NDEG)
        CHARC2=CHARC2/DFLOAT(NDEG)
        CHARINV=CHARINV/DFLOAT(NDEG)
        IF(DABS(CHARINV-1.0D0).LT.TOL)IRREP(3:3)='g'
        IF(DABS(CHARINV+1.0D0).LT.TOL)IRREP(3:3)='u'
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARCN-1.D0).LT.TOL)IRREP(1:1)='A'
         IF(DABS(CHARCN+1.D0).LT.TOL)IRREP(1:1)='B'
         IF(DABS(CHARC2-1.D0).LT.TOL)IRREP(2:2)='1'
         IF(DABS(CHARC2+1.D0).LT.TOL)IRREP(2:2)='2'
        ELSEIF(NDEG.EQ.2)THEN
         DO 402 I=1,ITOP
          IRREP(1:1)='E'
          Z=DCOS(ANGLE(I))
          IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(2:2),'(I1)')I
402      CONTINUE
        ELSE
         IRREP(2:2)='E'
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C S(2n) GROUPS.  SEPARATE CODE FOR n ODD AND EVEN
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'S')THEN
       IORDER=ATOI(PTGRP(2:ICLST))
       ITOP=(IORDER-1)/2
       DO 500 I=1,ITOP
        ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
500     CONTINUE
       IF(MOD(IORDER,2).NE.0)THEN
        CHARSN=CHRVEC(1)
        CHARINV=CHRVEC(IORDER/2)
        CHARSN=CHARSN/DFLOAT(NDEG)
        CHARINV=CHARINV/DFLOAT(NDEG)
        IF(DABS(CHARINV-1.0D0).LT.TOL)IRREP(4:4)='g'
        IF(DABS(CHARINV+1.0D0).LT.TOL)IRREP(4:4)='u'
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARSN-1.D0).LT.TOL)IRREP(3:3)='A'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(2:2)='E'
          DO 501 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(3:3),'(I1)')I
501       CONTINUE
         ELSE
          IRREP(3:3)='E'
         ENDIF
        ENDIF
C
       ELSE 
C
        CHARSN=CHRVEC(1)
        CHARSN=CHARSN/DFLOAT(NDEG)
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHARSN-1.D0).LT.TOL)IRREP(3:3)='A'
         IF(DABS(CHARSN+1.D0).LT.TOL)IRREP(3:3)='B'
        ELSEIF(NDEG.EQ.2)THEN
         IF(ITOP.NE.1)THEN
          IRREP(3:3)='E'
          DO 502 I=1,ITOP
           Z=DCOS(ANGLE(I))
           IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(4:4),'(I1)')I
502       CONTINUE
         ELSE
          IRREP(4:4)='E'
         ENDIF
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C C(n) GROUPS.
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'C'.AND.PTGRP(ICLST:ICLST).NE.'h'.AND.
     &       PTGRP(ICLST:ICLST).NE.'v')THEN
       IORDER=ATOI(PTGRP(2:ICLST))
       ITOP=IORDER/2
       DO 600 I=1,ITOP
        ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
600    CONTINUE
       CHARCN=CHRVEC(1)
       CHARCN=CHARCN/DFLOAT(NDEG)
       IF(NDEG.EQ.1)THEN
        IF(DABS(CHARCN-1.D0).LT.TOL)IRREP(4:4)='A'
        IF(DABS(CHARCN+1.D0).LT.TOL)IRREP(4:4)='B'
       ELSEIF(NDEG.EQ.2)THEN
        IF(ITOP.NE.1)THEN
         IRREP(3:3)='E'
         DO 601 I=1,ITOP
          Z=DCOS(ANGLE(I))
          IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(4:4),'(I1)')I
601      CONTINUE
        ELSE
         IRREP(4:4)='E'
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C D(n) GROUPS.  DOES NOT WORK FOR D2.
C**********************************************************************
C
      ELSEIF(PTGRP(1:1).EQ.'D'.AND.PTGRP(ICLST:ICLST).NE.'h'.AND.
     &       PTGRP(ICLST:ICLST).NE.'d')THEN
       IORDER=ATOI(PTGRP(2:ICLST))
       ITOP=IORDER/2
       DO 700 I=1,ITOP
        ANGLE(I)=DFLOAT(I)*TWOPI/DFLOAT(IORDER)
700    CONTINUE
       CHARCN=CHRVEC(1)
       CHARC2=CHRVEC(IORDER)
       CHARCN=CHARCN/DFLOAT(NDEG)
       CHARC2=CHARC2/DFLOAT(NDEG)
       IF(NDEG.EQ.1)THEN
        IF(DABS(CHARCN-1.D0).LT.TOL)IRREP(3:3)='A'
        IF(DABS(CHARCN+1.D0).LT.TOL)IRREP(3:3)='B'
        IF(DABS(CHARC2-1.D0).LT.TOL)IRREP(4:4)='1'
        IF(DABS(CHARC2+1.D0).LT.TOL)IRREP(4:4)='2'
       ELSEIF(NDEG.EQ.2)THEN
        IF(ITOP.NE.1)THEN
         IRREP(3:3)='E'
         DO 701 I=1,ITOP
          Z=DCOS(ANGLE(I))
          IF(DABS(CHARCN-Z).LT.TOL)WRITE(IRREP(4:4),'(I1)')I
701      CONTINUE
        ELSE
         IRREP(4:4)='E'
        ENDIF
       ENDIF
       RETURN
C
C**********************************************************************
C  THE FOLLOWING CODE PERTAINS TO THE CUBIC GROUPS T, Td, O, Th, Oh,
C  I AND Ih.  STANDARD ORIENTATIONS FOR THESE GROUPS ARE DESCRIBED IN
C  THE ACES II PROGRAM MANUAL 
C**********************************************************************
C
      ELSEIF(PTGRP(1:3).EQ.'T d')THEN
       IF(NDEG.EQ.2)THEN
        IRREP='   E'
        RETURN
       ELSE
        CHAR=CHRVEC(1)
        IF(NDEG.EQ.1)THEN
         IF(DABS(CHAR+1.D0).LT.TOL)IRREP='  A2'
         IF(DABS(CHAR-1.D0).LT.TOL)IRREP='  A1'
        ELSE
         IF(DABS(CHAR+1.D0).LT.TOL)IRREP='  T2'
         IF(DABS(CHAR-1.D0).LT.TOL)IRREP='  T1'
        ENDIF
       ENDIF
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'O h')THEN
       CHAR=CHRVEC(47)
       Z=CHAR/DFLOAT(NDEG)
       IF(DABS(Z+1.D0).LT.TOL)IRREP(4:4)='u'
       IF(DABS(Z-1.D0).LT.TOL)IRREP(4:4)='g'
       IF(NDEG.EQ.2)THEN
        IRREP(2:3)=' E'
        RETURN
       ENDIF
       Z1=CHRVEC(1)
       IF(DABS(Z1+1.D0).LT.TOL)IRREP(3:3)='2'
       IF(DABS(Z1-1.D0).LT.TOL)IRREP(3:3)='1'
       IF(NDEG.EQ.3)IRREP(2:2)='T'
       IF(NDEG.EQ.1)IRREP(2:2)='A'
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'O  ')THEN
       CHAR=CHRVEC(1)
       IF(DABS(CHAR-1.D0).LT.TOL)IRREP(4:4)='1'
       IF(DABS(CHAR+1.D0).LT.TOL)IRREP(4:4)='2'
       IF(NDEG.EQ.1)IRREP(3:3)='A'
       IF(NDEG.EQ.2)IRREP='   E'
       IF(NDEG.EQ.3)IRREP(3:3)='T'
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'T  ')THEN
       IF(NDEG.EQ.1)IRREP='   A'
       IF(NDEG.EQ.2)IRREP='   E'
       IF(NDEG.EQ.3)IRREP='   T'
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'T h')THEN
       CHAR=CHRVEC(23)
       CHAR=CHAR/DFLOAT(NDEG)
       IF(DABS(CHAR+1.D0).LT.TOL)IRREP(4:4)='u'
       IF(DABS(CHAR-1.D0).LT.TOL)IRREP(4:4)='g'
       IF(NDEG.EQ.1)IRREP(1:3)='  A'
       IF(NDEG.EQ.2)IRREP(1:3)='  E'
       IF(NDEG.EQ.3)IRREP(1:3)='  T'
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'I h')THEN
       CHAR=CHRVEC(119)
       CHAR=CHAR/DFLOAT(NDEG)
       IF(DABS(CHAR+1.D0).LT.TOL)IRREP(4:4)='u'
       IF(DABS(CHAR-1.D0).LT.TOL)IRREP(4:4)='g'
       IF(NDEG.EQ.1)IRREP(3:3)='A'
       IF(NDEG.EQ.4)IRREP(3:3)='G'
       IF(NDEG.EQ.5)IRREP(3:3)='H'
       IF(NDEG.EQ.3)THEN
        ETAP=0.5D0*(1.D0+SQRT(5.D0))
        ETAM=0.5D0*(1.D0-SQRT(5.D0))
        CHAR=CHRVEC(1)
        IF(DABS(CHAR-ETAP).LT.TOL)IRREP(2:3)='T1'
        IF(DABS(CHAR-ETAM).LT.TOL)IRREP(2:3)='T2'
       ENDIF
       RETURN
      ELSEIF(PTGRP(1:3).EQ.'I  ')THEN
       IF(NDEG.EQ.1)IRREP(4:4)='A'
       IF(NDEG.EQ.2)IRREP(4:4)='E'
       IF(NDEG.EQ.4)IRREP(4:4)='G'
       IF(NDEG.EQ.5)IRREP(4:4)='H'
       IF(NDEG.EQ.3)THEN
        ETAP=0.5D0*(1.D0+SQRT(5.D0))
        ETAM=0.5D0*(1.D0-SQRT(5.D0))
        CHAR=CHRVEC(1)
        IF(DABS(CHAR-ETAP).LT.TOL)IRREP(3:4)='T1'
        IF(DABS(CHAR-ETAM).LT.TOL)IRREP(3:4)='T2'
       ENDIF
       RETURN
      ELSE
C
C ERROR HANDLING IF PTGRP IS NOT RECOGNIZED IN CODE ABOVE
C
       WRITE(6,2500)PTGRP
2500   FORMAT(T3,'@DETIRR-W, Point group ',A,' unknown.')
      ENDIF
C
      RETURN
      END
