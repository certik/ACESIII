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
      SUBROUTINE SYMUNQ(NATOM,IORDGP,IATNUM,SYMOPS,IPTR,
     &IATLST,MEMBER,TRNMAT,ORBPOP,IORBIT,PTGRP)
C
C ARRAYS:
C
C  MEMBER  -  GIVES ATOM LIST GROUPED ACCORDING TO ORBITS.
C  TRNMAT(I)- GIVES INDEX OF SYMOP WHICH MAPS THE REFERENCE ATOM FOR
C              THE ORBIT TO WHICH I BELONGS INTO I.
C  ORBPOP(I)- THE POPULATION OF ORBIT(I).
C
C
C THIS ROUTINE GROUPS THE ATOMS INTO SYMMETRY EQUIVALENT SETS.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "io_units.par"
      COMMON /LOCAL/ STSYM(MXATMS)
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
 
 
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      CHARACTER*4 PTGRP
      CHARACTER*8 STSYM,SITGRP
      INTEGER TRNMAT(NATOM),CHKINT,IATNUM(NATOM),
     &MEMBER(NATOM),ORBPOP(NATOM)
      LOGICAL FULGRP
      DIMENSION SYMOPS(9*IORDGP),IPTR(NATOM,IORDGP),IATLST(NATOM)
      IP(I)=1+9*(I-1)
      FULGRP=.FALSE.
      IF(FPGRP(2:2).EQ.'X'.AND.PTGRP(2:2).EQ.'8'.AND.FPGRP(1:1).EQ.
     &   PTGRP(1:1).OR.PTGRP.EQ.FPGRP)FULGRP=.TRUE.
      IORBIT=0
c      CALL IZERO(IATLST,NATOM)
c      CALL IZERO(TRNMAT,NATOM)
      do i = 1, natom
         iatlst(i) = 0
         trnmat(i) = 0
      enddo
C
C SET ALL DUMMY ATOM POSITIONS TO 999.
C
      DO 5 I=1,NATOM
       IF(IATNUM(I).EQ.0)IATLST(I)=999
5     CONTINUE
C
C FIND FIRST ZERO IN IATLST AND USE THIS AS THE POSITION FOR THE
C  REFERENCE ATOM OF ORBIT IORBIT.
C
      ICOUNT=0
 1    IREF=CHKINT(IATLST,NATOM,0)
      JCOUNT=0
      IF(IREF.EQ.0)GOTO 999
      IORBIT=IORBIT+1
C
C CREATE TRNMAT AND MEMBER LISTS FOR THIS PARTICULAR ORBIT.
C
      CALL FRMLST(NATOM,IREF,IORBIT,IORDGP,TRNMAT,IATLST,IPTR,
     &                  MEMBER,ICOUNT,JCOUNT)
      ORBPOP(IORBIT)=JCOUNT
      GOTO 1
C
C NOW WRITE OUT TRANSFORMATION MATRICES, POINTER LISTS AND TRNMAT
C  LIST TO DISK FOR LATER PICKUP.  PUT TRNMAT AND MEMBER ON A
C  SEPARATE RECORD SINCE THEY COULD BE CHANGED IN PIKXYZ IF THE
C  REFERENCE ATOM FOR A PARTICULAR ORBIT IS NOT THE FIRST MEMBER
C  IN Z-MATRIX ORDERING.
C
999   DO 153 I=1,IORDGP
C       CALL MTRAN2(SYMOPS(IP(I)),3)
153   CONTINUE
C
C WRITE THIS STUFF ONLY IF THE CALCULATION IS A FINITE DIFFERENCE
C  CALC. AND FDARC IS NOT THERE.
C
      IF(IVIB.EQ.2.AND.FULGRP)THEN
       OPEN(UNIT=15,FILE='FDARC',FORM='UNFORMATTED',
     &     STATUS='UNKNOWN')
       REWIND(15)
       WRITE(15)NATOM,IORDGP,IORBIT,ORBPOP,SYMOPS,IPTR
       WRITE(15)MEMBER,TRNMAT
       CLOSE(UNIT=15,STATUS='KEEP')
      ENDIF
C
C WRITE OUT THIS INFORMATION FOR FINDIF ONLY.
C
      IF(IPRNT.GE.5.OR.(IVIB.EQ.2.AND.FULGRP))
     &THEN
       WRITE(LUOUT,8080)IORBIT,PTGRP
8080   FORMAT(T3,'@SYMUNQ-I, There are ',I3,' orbits in the ',
     &          A,' point group.')
c       CALL IZERO(TRNMAT,NATOM)
       do i = 1, natom
          trnmat(i) = 0
       enddo

       WRITE(LUOUT,9000)
       WRITE(LUOUT,5002)
       WRITE(LUOUT,9000)
5002   FORMAT(T5,'Set',T11,' Site Group',T45,'Members of set')
      ENDIF
      ICOUNT=0
      IOFF=0
      DO 30 I=1,IORBIT
       NORBIT=ORBPOP(I)
 
C
C GET SITE GROUP FOR ORBIT
C
       STSYM(I)=SITGRP(PTGRP,IORDGP,ORBPOP(I))//'    '
       IF(STSYM(I)(1:3).EQ.'XXX')THEN
C
C THIS IS EITHER C2h OR Dnd, FOR WHICH THE SITE GROUP MUST
C   BE EITHER Cs OR C2.  FIGURE THIS OUT NOW.
C
        CHRMIN=10.D0
        DO 55 IOP=1,IORDGP
         IREFER=MEMBER(IOFF+1)
         IF(IPTR(IREFER,IOP).EQ.IREFER)THEN
          CHRMIN=MIN(CHRMIN,XTRACE(SYMOPS(IP(IOP)),3))
         ENDIF
55      CONTINUE
        STSYM(I)(1:4)='C2  '
        IF(CHRMIN+0.9D0.GT.0.D0)STSYM(I)(1:4)='C s '
       ENDIF
C
C GO ON.
C
       IF((IVIB.EQ.2.AND.FULGRP).OR.IPRNT.GT.5)
     & WRITE(LUOUT,5001)I,STSYM(I)(1:3),(MEMBER(IOFF+K),K=1,NORBIT)
5001   FORMAT(T5,I3,T14,A,(T30,15I3,/))
       IOFF=IOFF+ORBPOP(I)
 30   CONTINUE
      IF(IPRNT.GT.5.OR.(IVIB.EQ.2.AND.FULGRP))
     &WRITE(LUOUT,9000)
9000  FORMAT(72('-'))
      RETURN
      END
