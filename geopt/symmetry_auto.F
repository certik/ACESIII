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
#include "flags.h" 
      SUBROUTINE SYMMETRY_AUTO(SCRATCH, QTMP, NEWQ, IT, IDEGEN, 
     &                         ORIEN2, NOSILENT)
C
C     SUBROUTINE DOES SYMMETRY ANALYSIS AND PUTS MOLECULE IN AN
C     ORIENTATION THAT ACES2 CAN DEAL WITH
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "fnamelen.par"
#include "coord.com"

      DOUBLE PRECISION IT(3,3),CM(3),IV(3,3),NEWQ(NX),dtmp
      DOUBLE PRECISION RM(3,3),QTMP(NX),TATB(3)
      DOUBLE PRECISION SCRATCH(3*NX),MOLWT,ORIEN2(9)
      LOGICAL ALLDONE,ITWOAX, TSSEARCH, NOSILENT 

cSB START
      logical olddone,tetrah
cSB END

      LOGICAL OPTRES
      CHARACTER*3 SYMSTR(6)
      CHARACTER*4 some_string_func, JNKSTR, TMPGRP
      CHARACTER*2 XYZP(3)
      CHARACTER*1 XYZ(3),DOSTR(3),CHRTMP
      CHARACTER*(fnamelen) FNAME
      INTEGER ORDNEW,ORDOLD,IORGRP
      LOGICAL YESNO,SAXIS,SKIP,XYZIN,NWFINDIF
      DIMENSION NORD(3*MXATMS),IORBPOP(MXATMS)
      DIMENSION IORBREF(MXATMS)
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON /TOLERS/ SYMTOL,DEGTOL
      COMMON/RESTART2/OPTRES
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
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
     $     ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $     XYZTol
      COMMON /INPTYP/ XYZIN, NWFINDIF
C
#include "cbchar.com"
C
      CHARACTER*(5*MXATMS) ZSYMUNI
C
      INTEGER    LuArc
      PARAMETER (LuArc  = 77)
      CHARACTER*(*) ArcFil
      PARAMETER    (ArcFil = 'OPTARC ')
#include "io_units.par"
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA XYZ /'X','Y','Z'/
      DATA XYZP /'YZ','XZ','XY'/
      DATA ONE   / 1.0D+00/
      DATA ONEM  /-1.0D+00/
      DATA ZILCH / 0.0D+00/
      DATA IONE /1/
      IND(I,J) = 3*(I-1)+J
c
c
      CALL GETCREC(-1, 'JOBARC', "PTGP    ", 4, FPGRP) 
      TMPGRP=FPGRP
      FPGRP='    '
      DTOR=DACOS(ONEM)/180.D0
C
C INITIALIZE MILLIONS OF VARIABLES THAT CFT77 SAYS AREN'T INITIALIZED.
C
      OLDDONE=.FALSE.
      IHIGH=0
      IROT=0
      IREF=0
      IINV=0
      IDGRP=0
      ILINEAR=0
      ICOMP=0
      ICOUNT=0
      IDONE=0
      ISINV=0
      ICOMPQ=0
      IIAX=0
      ORDIS=0
      IERR=0
      ISET=0
      NUMB=0
      ICOMP2=0
      IHIGHX=0
      ISAXIS=0
      JCOMPX=0
      JCOMP=0
      ICOMPX=0
      ISKIP=0
      DOSTR(1)='X'
      DOSTR(2)='Y'
      DOSTR(3)='Z'
      TSSEARCH = (INR .EQ. 4 .OR. INR .EQ. 5 .OR. INR .EQ. 6)
c      CALL IZERO(IORBPOP,MXATMS)
c      CALL IZERO(IORBREF,MXATMS)
      do i = 1, mxatms
         iorbpop(i) = 0
         iorbref(i) = 0
      enddo

      IF(IDEGEN.EQ.0)THEN
C IDEGEN is zero for Abelian point groups. We can transfer the
C logic to statement 9000 (only symmetry operations
C that we need to deal with is rotation, reflection, and inversion).
C It looks like that this goto 9000 can not be easily avoided.
        GOTO 9000
      ENDIF
C This 777 entry point is used when the degeneracy found by
C inertia matrix was accidental.
  777 CONTINUE
C
C****************************************************************
C
C CODE FOR POINT GROUPS WITH DOUBLY DEGENERATE REPRESENTATIONS.
C
C****************************************************************
      IF(IDEGEN.EQ.1)THEN
C
C     IF DOUBLY DEGENERATE, CHECK TO MAKE SURE THAT THE PRINCIPAL AXES
C     ARE PARALLEL TO X,Y AND Z.  IF NOT, ROTATE ORIENTATION TO THIS
C     POINT.
C
        FPGRP = '    '
        IF(IPRNT .GE. 3 .AND. NOSILENT)WRITE(LUOUT,8002)
        DO 1541 I=1,3
 1541   IF(DABS(IT(I,I)).LT.SYMTOL)ILINEAR=1
        IF(ILINEAR.EQ.1)THEN
          IF(IPRNT .GE. 3 .AND. NOSILENT)WRITE(LUOUT,8004)
C
C Linear molecules can be handled by rotation, reflection,
C and inversion, so transfer the logic to 9000
C
          GOTO 9000
        ENDIF
C
C IDENTIFY HIGHEST ORDER ROTATIONAL AXIS (CHECK THROUGH 24 )
C
        DO 10001 I = 2,24
          ZANG = 360.D0/DFLOAT(I)
          DO 10000 J = 3,3
            CALL DOSYOP('C',I,1,J,NATOMS,NEWQ,SCRATCH,IERR,0)
            CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMP,
     &                   SYMTOL)
            IF(ICOMP.EQ.0)THEN
              IF (XYZIN) CALL IPUTREC(20, "JOBARC", "SYMEQUIV", 
     &                               NATOMS, NORD)
              IHIGH  = I
              IHIGHX = J
              IF(IPRNT.GE.4 .AND. NOSILENT)WRITE(LUOUT,3150)XYZ(J),I
 3150         FORMAT(T3,'@SYMMETRY-I, ',A1,' is a C[',I2,'] axis.')
            ENDIF
10000     CONTINUE
10001   continue
C
        IF (XYZIN) CALL IPUTREC(20, "JOBARC", "HIGHSTAX", 1, IHIGH) 
C
        IF(IPRNT .GE. 4 .AND. NOSILENT)WRITE(LUOUT,800)IHIGH,
     &                                 XYZ(IHIGHX)
  800   FORMAT(T3,'@SYMMETRY-I, The highest order rotational ',
     &       'axis is C[',i2,'] about ',A1,'.')
        IF(IHIGH.EQ.0)THEN
          IF (IPRNT .GE. 4 .AND. NOSILENT)  WRITE(6,9301)
 9301     FORMAT(T3,'@SYMMETRY-W, An accidental degeneracy is present.',
     &         ' Analysis might be wrong.')
C
          CALL PERPOP2('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,SYMTOL,
     &         IFOUND)
          IF(IFOUND.NE.1)THEN
            CALL PERPOP2('P',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,SYMTOL,
     &           IFOUND)
          ENDIF
          IDEGEN=0
          GOTO 777
        ENDIF
C
C CHECK FOR S(2n) AXIS
C
        CALL DOSYOP ('S',2*IHIGH,1,IHIGHX,NATOMS,NEWQ,SCRATCH,IERR,
     &       IMODE)
        CALL COMPARE(SCRATCH,NEWQ,ATMASS,NORD,NATOMS,ICOMP,SYMTOL)
        SAXIS=(ICOMP.EQ.0)
        IF (SAXIS) CALL IPUTREC(20, "JOBARC", "SYMEQUIV", NATOMS, NORD)
        IF(IPRNT.GT.10.AND.SAXIS .AND. NOSILENT)WRITE(6,9503)
 9503   FORMAT(T3,'@SYMMETRY-I, S(2n) axis found.')
C
C NOW CHECK FOR PERPENDICULAR C2 AXES.  THIS SEPARATES D AND C
C GROUPS.  THE LATTER ARE HANDLED FIRST.
C
        CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,SYMTOL,IFOUND)
        IF(IFOUND.EQ.1)THEN
          IF(IPRNT.GT.10 .AND. NOSILENT)WRITE(6,9501)
 9501     FORMAT(T3,'@SYMMETRY-I, Perpendicular C2 axis found.')
C
          FPGRP='DN  '
C
C CHECK FOR HORIZONTAL PLANE NOW
C
          CALL DOSYOP('P',1,1,3,NATOMS,NEWQ,SCRATCH,IERR,IMODE)
          CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMPX,SYMTOL)
C
          IF(ICOMPX.EQ.0)THEN
          IF (XYZIN) CALL IPUTREC(20, "JOBARC", "SAMEPLNE", NATOMS,
     &                           NORD)
            FPGRP='DNh '
            IF(IPRNT.GT.10 .AND. NOSILENT)WRITE(6,9502)
 9502       FORMAT(T3,'@SYMMETRY-I, Horizontal plane found.')
          ENDIF
          IF(SAXIS)FPGRP='DNd '
          IF (SAXIS .AND. XYZIN) THEN
             CALL DOSYOP('C', 8, 1, 3, NATOMS, NEWQ, SCRATCH, IERR, 0)
             CALL DCOPY(3*NATOMS, NEWQ, 1, SCRATCH(3*NATOMS+1), 1)
             CALL REFLECT(SCRATCH,NEWQ,NATOMS,1)
             CALL COMPARE(NEWQ, SCRATCH, ATMASS, NORD, NATOMS, ICOMP,
     &                    SYMTOL)
             CALL DCOPY(3*NATOMS, SCRATCH(3*NATOMS+1), 1, NEWQ, 1)
             CALL IPUTREC(20, "JOBARC", "SAMEPLNE", NATOMS,
     &                   NORD)
          ENDIF
        ENDIF
C
        IF(FPGRP.EQ.'    ')THEN
C
C THIS IS NOT A D GROUP, SO WE KNOW THAT IT IS EITHER CN, CNv, CNh
C OR MAYBE SN.  PROCEED WITH LOGIC.
C
C FIRST CHECK FOR HORIZONTAL PLANE - THIS ESTABLISHES THE POINT GROUP
C AS CNh
C
          CALL DOSYOP ('P',1,1,3,NATOMS,NEWQ,SCRATCH,IERR,IMODE)
          CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMP,SYMTOL)
          IF(ICOMP.EQ.0)THEN
            FPGRP='CNh '
          IF (XYZIN) CALL IPUTREC(20, "JOBARC", "SAMEPLNE", NATOMS,
     &                           NORD)
            IF(IPRNT.GT.10 .AND. NOSILENT)WRITE(6,9502)
            GOTO 9000
          ENDIF
C
C CNv OR CN ?  THIS CAN BE DETERMINED BY PRESENCE OR ABSENCE OF
C PLANES CONTAINING THE SYMMETRY AXIS.  CHECK FOR THESE NOW.
C
          CALL PERPOP('P',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,SYMTOL,IFOUND)
          IF(IFOUND.EQ.1)FPGRP='CNv '
        ENDIF
        IF(FPGRP.EQ.'    ')THEN
          FPGRP='CN  '
          IF(SAXIS)THEN
            FPGRP='SN  '
            IHIGH=2*IHIGH
          ENDIF
        ENDIF
      ENDIF
C****************************************************************
C
C CODE FOR CUBIC POINT GROUPS.  BY FAR THE MOST COMPLICATED PART.
C
C****************************************************************
      IF(IDEGEN.EQ.2)THEN
        ALLDONE=.FALSE.

cSB the following is no longer necessary, it is left only so that
c   the cartesian coordinates are exactly identical to the old
c   version of joda
        itwoax =.false.
        tetrah =.false.
cSB

        IF(IPRNT .GE. 3 .AND. NOSILENT)WRITE(LUOUT,8003)
C
C CHECK FOR INVERSION SYMMETRY.
C
        DO 14002 I = 1,NATOMS*3
14002   SCRATCH(I) = -NEWQ(I)
        CALL COMPARE(SCRATCH,NEWQ,ATMASS,NORD,NATOMS,ISINV,SYMTOL)
C
C NOW PICK A REFERENCE ATOM IN THE MOLECULE **WHICH DOES NOT LIE AT
C THE ORIGIN** AND THEN LOOP OVER ALL OTHER ATOMS HAVING THE SAME
C ATOMIC NUMBER AND COMPARABLE DISTANCE FROM THE COM.  CHECK FOR
C AXES OF ORDER 3, 4 AND 5.
C
        IREF=-1
        IOFF=1
        DO 5000 IATOM=1,NATOMS
          X=XDNRM2(3,NEWQ(IOFF),1)*ATMASS(IATOM)
          IF(X.GT.SYMTOL.AND.IREF.EQ.-1)THEN
            IREF=IATOM
            XREF=X
            ZREF=X/ATMASS(IATOM)
            IBOT=1+(IREF-1)*3
          ENDIF
          IOFF=IOFF+3
 5000   CONTINUE
C
C  HANDLE ATOMIC CALCULATION LOGIC
C
        IF(IREF.EQ.-1)THEN
          FPGRP='I h'
          GOTO 9000
        ENDIF
        DO 5001 IATOM=1,NATOMS
          IBOT2=1+(IATOM-1)*3
          X=XDNRM2(3,NEWQ(IBOT2),1)*ATMASS(IATOM)
          IF(ABS(X-XREF).LT.SYMTOL)THEN
            IF(IPRNT.GT.100)WRITE(LUOUT,50000)IREF,IATOM
50000       FORMAT(T3,'@SYMMETRY-I, Checking atoms ',I3,' and ',I3,'.')
C
C FIRST ROTATE MOLECULE SO THAT THE INTERATOMIC DISTANCE VECTOR
C   IS PARALLEL TO THE Z-AXIS, WITH BISECTOR ALONG X-AXIS.
C   THIS MEANS THAT THE TWO ATOMS LIE IN THE XZ PLANE.
C   STRUCTURE NOW IN SCRATCH(1).  HERE ONE HAS TO ALLOW FOR THE
C   POSSIBILITY THAT ATOMS I AND J ARE 180 DEGREES APART; IF THIS
C   IS THE CASE, JUST GO TO BOTTOM OF LOOP SINCE THESE DO NOT
C   DEFINE ANY AXIS UNIQUELY.
C
            CALL VEC(NEWQ(IBOT),NEWQ(IBOT2),SCRATCH(1),0)
            dtmp = xdot(3,SCRATCH(1),1,SCRATCH(1),1)
            DIST = DSQRT(dtmp)
            CALL VADD(SCRATCH(1),NEWQ(IBOT),NEWQ(IBOT2),NX,ONE)
C
C DEAL WITH THE LINEAR PROBLEM RIGHT NOW.
C
            BILEN = xdot(3,SCRATCH(1),1,SCRATCH(1),1)
            IF(BILEN.LT.1.D-12)GOTO 5001
            CALL SIAZ(SCRATCH(1),RM,1)
            CALL ZERO(SCRATCH,NX*3)
C           CALL MATMULV(SCRATCH(NX+1),NEWQ,RM,NATOMS,3,3)
            CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,3,ZILCH,
     &                 SCRATCH(NX+1),3)
C
C NOW YOU HAVE BISECTOR ALONG X.  ROTATE ABOUT X TO BRING THE TWO
C  ATOMS INTO POSITION PARALLEL TO Z!
C
            dtmp = xdot(2,SCRATCH(NX+IBOT+1),1,SCRATCH(NX+IBOT+1),1)
            DIP  = DSQRT(dtmp)
            ARGU=SCRATCH(NX+IBOT+2)/DIP
            ANGLE2=-DACOSX(ARGU,1.D-10)/DTOR
            IF(SCRATCH(NX+IBOT+1).GT.0.D0)ANGLE2=-ANGLE2
            CALL ROTM(1,ANGLE2,1,RM)
C           CALL MATMULV(SCRATCH,SCRATCH(NX+1),RM,NATOMS,3,3)
            CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,SCRATCH(NX+1),3,
     &                 ZILCH,SCRATCH,3)
C
C NOW FIND NEW VECTOR WHICH BISECTS THE TWO ATOMS IN QUESTION -
C   IT HAD BEST BE X!
C
            CALL VADD(TATB,SCRATCH(IBOT),SCRATCH(IBOT2),
     &           3,ONE)
C
C 1. TEST IF THE TWO ATOMS ARE CONNECTED BY SYMMETRY AXES, LOOPING
C     OVER POSSIBILITIES (2,3,4 AND 5).
C
            DO 5002 IAXORD=5,2,-1
              ANGIAX=360.D0/DFLOAT(IAXORD)
              ANGMG=ANGMAG(ZREF,DIST,IAXORD,IERR)
              IF(IERR.EQ.0.AND..NOT.ALLDONE)THEN
                IF(IPRNT.GE.3 .AND. NOSILENT)WRITE(LUOUT,50001)ANGMG
50001           FORMAT(T3,' ANGMAG angle is ',f8.4)
                IF(IPRNT.GE.13 .AND. NOSILENT)WRITE(LUOUT,50002)IERR,IAXORD
50002           FORMAT(T3,' IERR is ',i2,' for ',i1)
                IF(IPRNT.GE.200 .AND. NOSILENT)THEN
                  WRITE(LUOUT,*)' IN ROTATIONAL LOOP, LOOKING FOR ',
     &                 'ORDER ',IAXORD
                  WRITE(LUOUT,50003)IREF,IATOM
50003             FORMAT(' PLAYING WITH ATOMS ',I2,' AND ',I2)
                  WRITE(LUOUT,*)' IBOT AND IBOT2 ARE ',IBOT,IBOT2
                  WRITE(LUOUT,*)' ORBIT DISTANCE IS ',ORDIS
                  WRITE(LUOUT,*)' INTERATOMIC DISTANCE IS ',DIST
                  WRITE(LUOUT,*)' MAGIC ANGLE IS ',ANGMG
                  WRITE(LUOUT,*)' MAGIC VECTOR IS ',(TATB(JJ),JJ=1,3)
                ENDIF
                CALL ROTM(3,ANGMG,1,RM)
                IF(IPRNT.GE.200 .AND. NOSILENT)WRITE(LUOUT,80)
     &            (SCRATCH(JZ),JZ=1,NX)
C               CALL MATMULV(SCRATCH(NX+4),TATB,RM,1,3,3)
                CALL XGEMM('T','N',3,1,3,ONE,RM,3,TATB,3,ZILCH,
     &                     SCRATCH(NX+4),3)
                CALL SIAZ(SCRATCH(NX+4),RM,3)
                CALL ZERO(SCRATCH(NX+1),NATOMS*3*2)
C               CALL MATMULV(SCRATCH(NX+1),SCRATCH,RM,NATOMS,3,3)
                CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,SCRATCH,3,ZILCH,
     &                     SCRATCH(NX+1),3)
                CALL ROTM(3,ANGIAX,1,RM)
C               CALL MATMULV(SCRATCH(2*NX+1),SCRATCH(NX+1),
C    &               RM,NATOMS,3,3)
                CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,SCRATCH(NX+1),
     &                     3,ZILCH,SCRATCH(2*NX+1),3)
                call compare (scratch(nx+1),scratch(2*nx+1),atmass,
     &               nord,natoms,iiax,symtol)
cOLD                WRITE(6,*) "The first one"
cOLD                WRITE(6,*) (NORD(I), I=1, NATOMS) 
                IF(XYZIN) CALL IPUTREC(20, "JOBARC", "SYMEQUIV", NATOMS, 
     &                                NORD)
                IF(IIAX.EQ.0)THEN
cSB START
c The following assumes that we will run into an O or I rotation before
c running into at least 1 C2 and 1 C3 rotation.  An example where this
c is not the case is an Oh molecule defined as two superimposed tetrahedrons
c (one the inverse of the other)
cSSS                  IF(IAXORD.EQ.2.AND.TETRAH)THEN
cSSS                    ALLDONE=.TRUE.
cSSS                    FPGRP='T  '
cSSS                    IF(ISINV.EQ.0)FPGRP='T h'
cSSS                    CALL XDCOPY (3*NATOMS,SCRATCH(NX+1),1,NEWQ,1)
cSSS                    CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,
cSSS     &                   SYMTOL,IFOUND)
cSSS                    IF(IFOUND.EQ.0)THEN
cSSS                      WRITE(6,50004)FPGRP
cSSS                      CALL ERREX
cSSS                    ENDIF
cSSS                    CALL DOSYOP('S',4,1,3,NATOMS,NEWQ,SCRATCH,IERR,0)
cSSS                    CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMP,
cSSS     &                   SYMTOL)
cSSS                    IF(ICOMP.EQ.0.AND.ISINV.NE.0)THEN
cSSS                      FPGRP='T d'
cSSSC
cSSSC COMMENT OUT TWO LINES BELOW TO GET Td TO RUN AS D2.
cSSSC
cSSS                      CALL DOSYOP('C',8,1,3,NATOMS,NEWQ,SCRATCH,IERR,0)
cSSS                      CALL XDCOPY (NX,SCRATCH,1,NEWQ,1)
cSSS                    ENDIF
cSSS                  ELSEIF(IAXORD.EQ.2.AND..NOT.TETRAH)THEN
cSSS                    ITWOAX=.TRUE.
cSSS                    CALL XDCOPY(3*NATOMS,SCRATCH(NX+1),1,QTMP,1)
cSSS                  ELSEIF(IAXORD.EQ.3)THEN
cSSS                    TETRAH=.TRUE.
cSSS                    IF(ITWOAX)THEN
cSSS                      ALLDONE=.TRUE.
cSSS                      FPGRP='T  '
cSSS                      IF(ISINV.EQ.0)FPGRP='T h'
cSSS                      CALL XDCOPY (NX,QTMP,1,NEWQ,1)
cSSS                      CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,
cSSS     &                     SYMTOL,IFOUND)
cSSS                      IF(IFOUND.EQ.0)THEN
cSSS                        WRITE(6,50004)FPGRP
cSSS                        CALL ERREX
cSSS                      ENDIF
cSSS                      CALL DOSYOP('S',4,1,3,NATOMS,NEWQ,SCRATCH,IERR,0)
cSSS                      CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,
cSSS     &                     ICOMP,SYMTOL)
cSSS                      IF(ICOMP.EQ.0.AND.ISINV.NE.0)FPGRP='T d'
cSSS                    ENDIF
c We replace the above with the following which essentiall ignores all
c C2 and C3 rotations.  If, when we are done, we have not encountered a
c O or I rotation, we know it must be a T type molecule.
c
c All of the logic involving olddone, tetrah, and iaxord is optional.
c It is included only to make sure that the cartesian coordinates are
c exactly identical (as opposed to being off by a symmetry rotation).
                  if (iaxord.eq.2 .and. .not.olddone) then
                    call xdcopy(3*natoms,scratch(nx+1),1,qtmp,1)
                    if (tetrah) then
                      olddone=.true.
                    else
                      itwoax=.true.
                    endif
                  elseif (iaxord.eq.3 .and. .not.olddone) then
                    tetrah=.true.
                    if (itwoax) olddone=.true.
cSB END

                  ELSEIF(IAXORD.EQ.5)THEN
                    FPGRP='I  '
                    IF(ISINV.EQ.0)FPGRP='I h'
                    CALL XDCOPY (3*NATOMS,SCRATCH(NX+1),1,NEWQ,1)
C
C IF WE HAVE FOUND A FIVE-FOLD AXIS, FIRST REORIENT MOLECULE SO
C THAT C2 AXES LIE ALONG X, Y AND Z
C
                    CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,
     &                   SYMTOL,IFOUND)
                    IF(IFOUND.EQ.0)THEN
                      WRITE(6,50004)FPGRP
                      CALL ERREX
                    ENDIF
                    CALL ROTM(2,90.D0,1,RM)
C                   CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
                    CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,
     &                         3,ZILCH,SCRATCH,3)
                    CALL XDCOPY (3*NATOMS,SCRATCH,1,NEWQ,1)
                    CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,
     &                   SYMTOL,IFOUND)
                    IF(IFOUND.EQ.0)THEN
                      WRITE(6,50004)FPGRP
                      CALL ERREX
                    ENDIF
                    CALL ROTM(2,90.D0,1,RM)
C                   CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
                    CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,
     &                         3,ZILCH,SCRATCH,3)
                    CALL XDCOPY (3*NATOMS,SCRATCH,1,NEWQ,1)
C
C NOW MAKE SURE THAT C5 ORIENTATION IS CORRECT
C
                    SQ5=ONE/DSQRT(5.D0)
                    ARGX=DSQRT(0.5D0*(ONE-SQ5))
                    ANGL=DACOS(ARGX)
                    CALL ROTM(1,ANGL,0,RM)
C                   CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
                    CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,
     &                         3,ZILCH,SCRATCH,3)
                    CALL ROTM(3,72.D0,1,RM)
C                   CALL MATMULV(SCRATCH(NX+1),SCRATCH,RM,NATOMS,3,3)
                    CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,
     &                         3,ZILCH,SCRATCH,3)
                    call compare(scratch(nx+1),scratch,atmass,
     &                   nord,natoms,icomp,symtol)
                    IF(ICOMP.NE.0)THEN
                      CALL ROTM(3,90.D0,1,RM)
C                     CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
                      CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,
     &                         3,ZILCH,SCRATCH,3)
                      CALL XDCOPY(NX,SCRATCH,1,NEWQ,1)
                    ENDIF
                    ALLDONE=.TRUE.
                  ELSEIF(IAXORD.EQ.4)THEN
C
C IF WE HAVE FOUND A FOUR-FOLD AXIS, THEN WE ARE DONE.  THE MOLECULE
C BELONGS TO EITHER O OR Oh AND REORIENTATION IS COMPLETELY
C STRAIGHTFORWARD.
C
cSSS                    FPGRP='O  '
cSSS                    IF(ISINV.EQ.0)FPGRP='O h'
                    CALL XDCOPY (3*NATOMS,SCRATCH(NX+1),1,NEWQ,1)
                    CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,
     &                   SYMTOL,IFOUND)
                    IF(IFOUND.EQ.0)THEN
                      WRITE(6,50004)FPGRP
50004                 FORMAT(T3,'@SYMMETRY-F, Reorientation failure ',
     &                     'for point group ',A,'.')
                      CALL ERREX
                    ENDIF

cSB START
c At this point, we have either an O, Oh, or Th molecule.  In either
c case, the molecule is oriented in one of two ways (visualized with a
c cube).
c   A.  4 vertices in the XZ plane, 4 in the YZ plane (X and Y axis
c       point to edges of a cube, Z axis points to a face)
c   B.  45 degree rotation from this (all axis point to faces)
c In the B orientation, there is no way (using an X, Y, Z rotation or
c an XY, XZ, YZ reflection) to distinguish between Oh and Th.  In The
c A orientation, an XY, XZ, YZ reflection fails for Th.  Also, in the
c B orientation (but NOT the A orientation), there is a C4 rotation
c along the X axis.
c
c If the molecule is O or Oh, we want to rotate it to the B orientation
c in the end.  If it is a Th molecule, we need to exit without having
c modified the orientation.  In this case, the above XDCOPY and PERPOP
c may cause a problem (but it may not).

c Look for the C4(x) rotation (icomp=0 if B orientation))
                    CALL DOSYOP('C',4,1,1,NATOMS,NEWQ,SCRATCH,IERR,0)
                    CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMP,
     &                   SYMTOL)

c The B orientation.
                    if (icomp.eq.0) then

c      In the B orientation, an inversion center means either Oh or Th.
c      Rotate to the A orientation and check for a any reflection plane.
                      if (isinv.eq.0) then
                        call dosyop('C',8,1,3,natoms,newq,scratch,
     &                       ierr,0)
                        call dosyop('P',1,1,3,natoms,scratch,
     &                       scratch(nx+1),ierr,0)
                        call compare(scratch,scratch(nx+1),atmass,
     &                       nord,natoms,icomp,symtol)
cSSS                        WRITE(6,*) "The First one"
cSSS                        WRITE(6,*) (NORD(II), II =1, NATOMS) 
cSSS                        CALL IPUTREC(20, "JOBARC", "SYMEQUIV", NATOMS, 
cSSS     &                              NORD)
                        if (icomp.eq.0) then
                          fpgrp='O h'
                          alldone=.true.
                        endif

c      In the B orientation, no inversion center means O.
                      else
                        fpgrp='O  '
                        alldone=.true.
                      endif

c The A orientation.
                    else

c      In the A orientation, an inversion center means either Oh or Th.
c      We need to check for a reflection plane.  If found, rotate to
c      B orientation.
                      if (isinv.eq.0) then
                        call dosyop('P',1,1,3,natoms,newq,scratch(nx+1),
     &                       ierr,0)
                        call compare(newq,scratch(nx+1),atmass,nord,
     &                       natoms,icomp,symtol)
cSSS                        WRITE(6,*) "The Second One"
cSSS                        WRITE(6,*) (NORD(II), II =1, NATOMS) 
cSSS                        CALL IPUTREC(20, "JOBARC", "SYMEQUIV", NATOMS, 
cSSS     &                              NORD)
                        if (icomp.eq.0) then
                          call dosyop('C',8,1,3,natoms,newq,scratch,
     &                         ierr,0)
                          call xdcopy (nx,scratch,1,newq,1)
                          fpgrp='O h'
                          alldone=.true.
                        endif

c      In the A orientation, no inversion center means O.  Rotate it
c      to the B orientation.
                      else
                        call dosyop('C',8,1,3,natoms,newq,scratch,
     &                       ierr,0)
                        call xdcopy (nx,scratch,1,newq,1)
                        fpgrp='O  '
                        alldone=.true.
                      endif
                    endif

c The old logic was wrong
cSSS              CALL DOSYOP('C',4,1,1,NATOMS,NEWQ,SCRATCH,IERR,0)
cSSS              CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMP,
cSSS     &                     SYMTOL)
cSSS              IF(ICOMP.NE.0)THEN
cSSS               CALL DOSYOP('C',8,1,3,NATOMS,NEWQ,SCRATCH,IERR,0)
cSSS               CALL XDCOPY (NX,SCRATCH,1,NEWQ,1)
cSSS               ALLDONE=.TRUE.
cSSS              ENDIF
cSSSCJDW 12/22/94
cSSSC     This lets the Cr(CO)6 job run. CrC6 ends up with ICOMP=6 and so
cSSSC     ALLDONE is set to true. The C coordinates are different in Cr(CO)6
cSSSC     and CrC6. In the former all C atoms lie along Cartesian axes, while
cSSSC     in the latter four of the C atoms bisect axes. This is a somewhat
cSSSC     empirical fix, and a rigorous explanation would be appreciated.
cSSS              ALLDONE = .TRUE.
cSB END

                  ENDIF
                ENDIF
              ENDIF
 5002       CONTINUE
          ENDIF
 5001   CONTINUE
        IF(.NOT.ALLDONE)THEN

cSB START
c It used to be that if we did not get the ALLDONE flag, there was
c an error.  Now it just means that we are in a T-type group.

cSSS          WRITE(6,50005)
cSSS50005     FORMAT(T3,'@SYMMETRY-F, Cubic group not determined.')
cSSS          CALL ERREX

          FPGRP='T  '
          IF(ISINV.EQ.0)FPGRP='T h'
          CALL XDCOPY (NX,QTMP,1,NEWQ,1)
          CALL PERPOP('C',NEWQ,SCRATCH,ATMASS,NATOMS,NORD,
     &         SYMTOL,IFOUND)
          IF(IFOUND.EQ.0)THEN
            WRITE(6,50004)FPGRP
            CALL ERREX
          ENDIF
          CALL DOSYOP('S',4,1,3,NATOMS,NEWQ,SCRATCH,IERR,0)
          CALL COMPARE(NEWQ,SCRATCH,ATMASS,NORD,NATOMS,ICOMP,
     &         SYMTOL)
          IF(ICOMP.EQ.0.AND.ISINV.NE.0)THEN
            FPGRP='T d'
            IF (XYZIN) CALL IPUTREC(20, "JOBARC", "SYMEQUIV", NATOMS, 
     &                             NORD)
C
C COMMENT OUT TWO LINES BELOW TO GET Td TO RUN AS D2.
C
            CALL DOSYOP('C',8,1,3,NATOMS,NEWQ,SCRATCH,IERR,0)
            CALL XDCOPY (NX,SCRATCH,1,NEWQ,1)
          ENDIF
cSB END

        ENDIF
      ENDIF
C
C CUBIC POINT GROUP NOW DETERMINED.  ON WITH THE SHOW.
C
 80   FORMAT((4X,3(2X,F16.12)))
 8001 FORMAT(T3,'@SYMMETRY-I, The molecule belongs to an ',
     &     'Abelian group.')
 8002 FORMAT(T3,'@SYMMETRY-I, The molecule belongs to a point ',
     &     'group with doubly',
     &     ' degenerate representations.')
 8004 FORMAT(T3,'@SYMMETRY-I, The molecule is linear.')
 8003 FORMAT(T3,'@SYMMETRY-I, The molecule belongs to a cubic ',
     &     'point group.')
C
C This is the principal axis or computational orientation of the molecule
C
      IF(OPTRES) THEN
        CALL DGETREC(20,'JOBARC','COORD   ',NX,NEWQ)
      ENDIF
C
C This is the entry point for Abelain subgrop. All these 'GO TO 9000"
C could have been easily avoided.
 9000 CONTINUE
      IF(IPRNT.GE.0 .AND. NOSILENT) THEN
        WRITE(LUOUT,8202)
 8202   FORMAT(T3,' Principal axis orientation for molecule:')
        WRITE(LUOUT,80)(NEWQ(I),I=1,NX)
      ENDIF
C
C     CHECK SYMMETRY OPERATIONS BELONGING TO ABELIAN GROUPS
C
C  THE SYMSTR INFORMATION IS USED ONLY TO MAKE VMOL DECKS.  COMMENTED
C   OUT FOR NOW SINCE SYMEQV AND GRPOPS ASSUME "FUDGED" ORIENTATIONS.
C
      ANG = 180.D0
      IREF=0
      IOPS=1
      DO 873 I=1,6
        SYMSTR(I)='   '
  873 CONTINUE
      DO 31 I = 3,1,-1
        CALL REFLECT(NEWQ,SCRATCH,NATOMS,I)
        call compare(newq,scratch,atmass,
     &       nord,natoms,icomp,symtol)
        IF(ICOMP.EQ.0)THEN
C               SYMSTR(IOPS)(2:2)=XYZ(I)
          IF (XYZIN) THEN
              CALl IGETREC(0, "JOBARC", "SAMEPLNE", LENGTH, IJUNK)
              IF (LENGTH .LT. 0) CALL IPUTREC(20, "JOBARC", "SAMEPLNE",
     &                                   NATOMS, NORD)
          ENDIF
          IOPS=IOPS+1
          IF(IPRNT .GE. 3 .AND. NOSILENT)WRITE(LUOUT,77)I
   77     FORMAT(T3,' Reflection in plane ',i2,' is a valid ',
     $         'symmetry operation.')
          IREF = IREF + 2**I/2
        ENDIF
        CALL ROTM(I,ANG,1,RM)
        IF(IPRNT .GE. 5 .AND. NOSILENT)WRITE(LUOUT,80)
     $       ((RM(IX,JX),JX = 1,3),IX = 1,3)
C       CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
        CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,3,ZILCH,SCRATCH,3)
        call compare(newq,scratch,atmass,
     &       nord,natoms,icomp,symtol)
        IF(ICOMP.EQ.0)THEN
C              SYMSTR(IOPS)(2:3)=XYZP(I)
          IOPS=IOPS+1
          IF(IPRNT .GE. 3 .AND. NOSILENT)WRITE(LUOUT,78)I
   78     FORMAT(T3,' Rotation about ',i2,' is a valid symmetry ',
     $         'operation ')
          IROT = IROT + 2**I/2
        endif
   31 continue
      DO 102 I = 1,NATOMS*3
  102 SCRATCH(I) = -NEWQ(I)
      call compare(newq,scratch,atmass,
     &     nord,natoms,icomp,symtol)
      IF(ICOMP.EQ.0)THEN
C            SYMSTR(IOPS)='XYZ'
        IOPS=IOPS+1
        IF(IPRNT .GE. 3 .AND. NOSILENT)WRITE(LUOUT,81)
   81   FORMAT(T3,' The molecule possesses an inversion center. ')
        IF (XYZIN) CALL IPUTREC(20, "JOBARC", "SYMEQUIV",
     &                              NATOMS, NORD)
        IF(ILINEAR.EQ.1)FPGRP = 'DXh '
        IINV=1
      ELSE
        IF(ILINEAR.EQ.1)FPGRP = 'CXv '
cOLD        IF (XYZIN) CALL IPUTREC(20, "JOBARC", "SYMEQUIV",
cOLD     &                              NATOMS, NORD)    
        CALL IGETREC(0, "JOBARC", "SYMEQUIV", LENGTH, IJUNK)
        IF (XYZIN .AND. LENGTH .LT. 0) CALL IPUTREC(20, "JOBARC", 
     &                                 "SYMEQUIV", NATOMS, NORD)    
      ENDIF
C
C     DETERMINE LARGEST ABELIAN SUBGROUP OF MOLECULE
C
      IF(IROT.EQ.0.AND.IREF.EQ.0.AND.IINV.EQ.0)PGRP = 'C1  '
      IF(IROT.EQ.0.AND.IREF.EQ.0.AND.IINV.EQ.1)PGRP = 'C i '
      IF(IROT.EQ.0.AND.IREF.NE.0.AND.IINV.EQ.0)PGRP = 'C s '
      IF(IROT.NE.0.AND.IREF.EQ.0.AND.IINV.EQ.0)PGRP = 'C2  '
      IF(IROT.NE.0.AND.IREF.NE.0.AND.IINV.EQ.0)PGRP = 'C2v '
      IF(IROT.NE.0.AND.IREF.NE.0.AND.IINV.EQ.1)PGRP = 'C2h '
      IF(IROT.EQ.7.AND.IREF.EQ.0.AND.IINV.EQ.0)PGRP = 'D2  '
c         IF(IROT.EQ.7.AND.IREF.EQ.0.AND.IINV.EQ.0)PGRP = 'C2  '
      IF(IROT.EQ.7.AND.IREF.EQ.7.AND.IINV.EQ.1)PGRP = 'D2h '
      IF(IPRNT .GE. 4 .AND. NOSILENT)WRITE(LUOUT,27)IREF,IROT,IINV
   27 FORMAT(T3,'Symmetry bits: ',3(1x,i3))
CJDW 9/16/96. Comment out 3 lines below, as they will not be needed if
C             JFS modifications work (routine SETD2XYZ). Above C2 line
C             is how JFS used to handle D2, i.e. switch to C2 always
C             was done.
CJDW 4/ 1/97. Putting these three lines back in. D2 is causing too much
C             trouble in finite difference calculations at the moment.
c
CNO 5/11/94 Remove D2 from findif calculations and also from 
C Resonance Raman Calculations.
C This restriction should no longer be needed
C (see my comments in vmlgen in vmol2ja and vib1 and geopt in joda)
C Ajith Perera, 01/2006.
CSSS
CSSS         if(pgrp.eq.'D2  ') then
CSSS           if(iflags(54).eq.3.or.iflags(54).eq.2.or.
CSSS     &        iflags2(3) .ne. 0)    pgrp = 'C2  '
CSSS         endif
c
C IFLAGS(85) is non-zero when the subgroup key-word is turned on
C
      IF(IFLAGS(85).NE.0)THEN
        IF(IPRNT .GE. 4 .AND. NOSILENT) WRITE(6,1051)
 1051   FORMAT(T3,'@SYMMETRY-I, Reference coordinates used for ',
     &       'subgroup specification.')
        CALL DUMPCORD(NATOMS,NEWQ,IATNUM)
        IF(FPGRP.EQ.'    ')FPGRP=PGRP
        IF(ISYM .EQ. 1)    PGRP='C1  '
#ifdef _DEBUG_LVL0
        Print*, "PGRP, FPGRP:", FPGRP, PGRP
#endif
        CALL SUBGROUP(IFLAGS(85),PGRP,FPGRP,IROTATE)
C 
        IF(IROTATE.EQ.1)THEN
          CALL ROTM(3,45.0D0,1,RM)
C         CALL MATMULV(SCRATCH,NEWQ,RM,NATOMS,3,3)
          CALL XGEMM('T','N',3,NATOMS,3,ONE,RM,3,NEWQ,3,ZILCH,SCRATCH,3)
          CALL XDCOPY(3*NATOMS,SCRATCH,1,NEWQ,1)
        ENDIF
        IF(IFLAGS(86)+1.NE.1)THEN
          CHRTMP=DOSTR(1)
          DOSTR(1)=DOSTR(IFLAGS(86)+1)
          DOSTR(IFLAGS(86)+1)=CHRTMP
        ENDIF
C
C In order to do geo. optimization in a subgroup, Andrew Taube, 06/20
C
        ORDNEW=IORGRP(PGRP)                
        ORDOLD=IORGRP(FPGRP)              
        IF (ORDNEW.LE.ORDOLD) FPGRP=PGRP 
      ENDIF
      CALL ZERO(ORIENT,9)
C
C     ROTATE TO DEFAULT SYMMETRY FRAME IF NEEDED - COMPLETELY DONE BY
C     BRUTE FORCE.  NO CUTE ALGORITHM USED HERE.
C
C Reorientation of the molecule depending on the point group.
      BPGRP = PGRP
      CALL ZERO(SCRATCH,NATOMS*3)
      IF(PGRP.EQ.'C1  '.OR.PGRP.EQ.'C i '.OR.
     &     PGRP.EQ.'D2  '.OR.PGRP.EQ.'D2h ')THEN
c            IF(PGRP(1:2).NE.'D2')PGRP = 'C1  '
        CALL VADD(Q,NEWQ,SCRATCH(1),NATOMS*3,ONE)
        GOTO 75
      ENDIF

C
C     C2V
C
      IF(PGRP.EQ.'C2v'.OR.PGRP.EQ.'C2h'.OR.PGRP.EQ.'C2 ')THEN
C
C     IF X IS ROTATION AXIS - SWITCH TO Z - AND PUT LARGEST MOMENT OF
C     INERTIA AROUND X - EIGENVECTORS RETURNED FROM EIG ARE SORTED
C     HIGHEST TO LOWEST FOR D2H, WE TAKE STANDARD ORIENTATION SUCH
C     THAT IX>IY>IZ
C
        CALL SETXYZ(Q,NEWQ,ORIENT,DOSTR(1),IROT,IOK,NATOMS)
        IF(IOK.EQ.0.AND.IFLAGS(86).EQ.0)THEN
          CALL SETXYZ(Q,NEWQ,ORIENT,DOSTR(2),IROT,IOK,NATOMS)
          IF(IOK.EQ.0.AND.IFLAGS(86).EQ.0)THEN
            CALL SETXYZ(Q,NEWQ,ORIENT,DOSTR(3),IROT,IOK,NATOMS)
          ENDIF
        ENDIF
C
      ELSE
C
        CALL SETXYZ(Q,NEWQ,ORIENT,DOSTR(1),IREF,IOK,NATOMS)
        IF(IOK.EQ.0.AND.IFLAGS(86).EQ.0)THEN
          CALL SETXYZ(Q,NEWQ,ORIENT,DOSTR(2),IREF,IOK,NATOMS)
          IF(IOK.EQ.0.AND.IFLAGS(86).EQ.0)THEN
            CALL SETXYZ(Q,NEWQ,ORIENT,DOSTR(3),IREF,IOK,NATOMS)
          ENDIF
        ENDIF
C
      ENDIF
c            IF(Mod(IRot,2).EQ.1)THEN
c               DO 1001 J = 1,NATOMS
c                  Q(3*J) = NEWQ(3*J-2)
c                  Q(3*J-2) = NEWQ(3*J-1)
c 1001             Q(3*J-1) = NEWQ(3*J)
c               ORIENT(1,2)=ONE
c               ORIENT(2,3)=ONE
c               ORIENT(3,1)=ONE
c            ELSEIF(Mod(IRot/2,2).EQ.1)THEN
c               DO 1002 J = 1,NATOMS
c                  Q(3*J) = NEWQ(3*J-1)
c                  Q(3*J-1) = NEWQ(3*J-2)
c 1002             Q(3*J-2) = NEWQ(3*J)
c               ORIENT(3,2)=ONE
c               ORIENT(2,1)=ONE
c               ORIENT(1,3)=ONE
c            ELSEIF(Mod(IRot/4,2).EQ.1)THEN
c               CALL VADD(Q,NEWQ,SCRATCH(1),NATOMS*3,ONE)
c               ORIENT(1,1)=ONE
c               ORIENT(2,2)=ONE
c               ORIENT(3,3)=ONE
c            ENDIF
c         ELSE
c            IF(Mod(IRef,2).EQ.1)THEN
c               DO 2001 J = 1,NATOMS
c                  Q(3*J) = NEWQ(3*J-2)
c                  Q(3*J-2) = NEWQ(3*J-1)
c 2001             Q(3*J-1) = NEWQ(3*J)
c               ORIENT(3,1)=ONE
c               ORIENT(1,2)=ONE
c               ORIENT(2,3)=ONE
c            ELSEIF(Mod(IRef/2,2).EQ.1)THEN
c               DO 2002 J = 1,NATOMS
c                  Q(3*J) = NEWQ(3*J-1)
c                  Q(3*J-1) = NEWQ(3*J-2)
c 2002             Q(3*J-2) = NEWQ(3*J)
c               ORIENT(3,2)=ONE
c               ORIENT(2,1)=ONE
c               ORIENT(1,3)=ONE
c            ELSEIF(MOD(IRef/4,2).EQ.1)THEN
c               CALL VADD(Q,NEWQ,SCRATCH(1),NX,ONE)
c               ORIENT(1,1)=ONE
c               ORIENT(2,2)=ONE
c               ORIENT(3,3)=ONE
c            ENDIF
c         ENDIF
C
C NOW RESET ORIENT TO THE IDENTITY IF THE GROUP IS ABELIAN.
C
C This is important, for abelian groups, the ORIENT is the unit matrix
      IF(IDEGEN.EQ.0)THEN
        CALL ZERO(ORIENT,9)
        ORIENT(1,1)=ONE
        ORIENT(2,2)=ONE
        ORIENT(3,3)=ONE
      ENDIF
C
C FILL SYMSTR STRING AND WRITE OUT IF THIS IS A VMOL-BASED CALCULATION.
C
C
C     IF THIS IS A VMOL-DRIVEN CALCULATION, WRITE SYMMETRY INFORMATION
C      TO VMLSYM. ALSO FOR ARGOS INPUT
C
CJDW 9/20/96. Also for SEWARD.
CADY 5/06/04. Also for GAMESS.
C
   75 IF(INTTYP.EQ.1.OR.INTTYP.EQ.2.OR.INTTYP.EQ.0.OR.INTTYP.EQ.5
     &     .OR.INTTYP.EQ.4)THEN
        IF(PGRP.EQ.'C2v')THEN
          SYMSTR(1)(2:2)='X'
          SYMSTR(2)(2:2)='Y'
          IOPS=3
        ELSEIF(PGRP.EQ.'D2h')THEN
          SYMSTR(1)(2:2)='X'
          SYMSTR(2)(2:2)='Y'
          SYMSTR(3)(3:3)='Z'
          IOPS=4
        ELSEIF(PGRP.EQ.'C2h')THEN
          SYMSTR(1)(2:2)='Z'
          SYMSTR(2)(2:3)='XY'
          IOPS=3
        ELSEIF(PGRP.EQ.'C s')THEN
          SYMSTR(1)(2:2)='Z'
          IOPS=2
        ELSEIF(PGRP.EQ.'C i')THEN
          SYMSTR(1)     ='XYZ'
          IOPS=2
        ELSEIF(PGRP.EQ.'D2 ')THEN
CJDW 9/16/96. Modification from JFS.
          CALL SETD2XYZ(NATOMS,Q,ATMASS,SYMSTR(1),SYMSTR(2),IERROR)
          IOPS=3
          IF(IERROR.EQ.1)THEN
            PGRP='C2 '
            SYMSTR(1)(2:3)='XY'
            IOPS=2
          ENDIF
c           SYMSTR(1)(2:3)='XY'
c           SYMSTR(2)(2:3)='XZ'
c           IOPS=3
        ELSEIF(PGRP.EQ.'C2 ')THEN
          SYMSTR(1)(2:3)='XY'
          IOPS=2
        ELSEIF(PGRP.EQ.'C1 ')THEN
          IOPS=1
        ENDIF
        INQUIRE(FILE='VMLSYM',EXIST=YESNO,NUMBER=JUNK)
        IF(YESNO)THEN
          OPEN(UNIT=78,FILE='VMLSYM',STATUS='OLD',FORM='UNFORMATTED')
          CLOSE(UNIT=78,STATUS='DELETE')
        ENDIF
        OPEN(UNIT=78,FILE='VMLSYM',STATUS='NEW',FORM='UNFORMATTED')
        WRITE(78)MIN(3,IOPS-1),(SYMSTR(I),I=1,3)
        CLOSE(UNIT=78,STATUS='KEEP')
      ENDIF
C
C     DO NUMERIC-ASCII CONVERSION OF IHIGH.
C
      JnkStr=some_string_func(FPGRP,IHIGH)
      FPGrp = JnkStr
      IF(FPGRP.EQ.'    ')FPGRP=BPGRP
C
C FIGURE OUT WHAT THE COMPUTATIONAL POINT GROUP IS FROM OPTCTL PARAMETER
C  ISYM.  IF IT IS SET TO 0, THEN JUST USE PGRP AS IS.  IF =1, THEN USE
C  C1.  IF =3, USE FULL POINT GROUP{SYMM=NONE(0),OFF(1),ON(2),FULL(3)}.
C
      IF(ISYM.EQ.1)PGRP='C1  '
      IF(ISYM.EQ.3)PGRP=FPGRP
C The ncycle=0 refers to the first joda run during geo. optimizations
      IF ( ncycle .EQ. 0 ) THEN
        IF(NOSILENT)WRITE(LUOUT,788)
        IF(IDEGEN.GT.0)THEN
          IF(NOSILENT)WRITE(LUOUT,1881)FPGRP
        ENDIF
        IF(IDEGEN.EQ.0)THEN
          IF(NOSILENT) 
     &       WRITE(LUOUT,1881)some_string_func(BPGRP,IHigh)
          FPGRP=BPGRP
        ENDIF
 1881   FORMAT(T3,' The full molecular point group is ',a,'.')
        IF(NOSILENT) WRITE(LUOUT,177)BPGRP
        IF(NOSILENT) WRITE(LUOUT,712)PGRP
        IF(NOSILENT) WRITE(LUOUT,788)
  788   FORMAT(80('*'))
  712   FORMAT(T3,' The computational point group is ',A,'.')
  177   FORMAT(T3,' The largest Abelian subgroup of the full ',
     $       'molecular point group is ',a,'.')
      ELSE
C
C  MAKE SURE THAT SYMMETRY HAS NOT BEEN REDUCED FROM PREVIOUS STEP.  IF
C    IT HAS, AND THIS IS NOT A TS SEARCH, THEN ABORT.
C
        IF(FPGRP.NE.TMPGRP)THEN
          ORDNEW=IORGRP(FPGRP)
          ORDOLD=IORGRP(TMPGRP)
          WRITE(LUOUT,1883)FPGRP
 1883     FORMAT(T3,'@SYMMETRY-I, Point group has changed to ',A,'.')
C
C Originally test for TS searches were done by INR .NE. 2. That
C was dangerous. Now the TSSEARCH is  a logical variable that 
C is set to TRUE for TS searches (see above, INR=4, 5, and 6). 
C Ajith Perera 08/2001. 
C 
C A bug fix, we must allow only the TS searches to proceed when
C there is a symmetry lowering. It is amazing that it took us
C 15 years to get this right! Ajith Perera 05/2005
C
        IF (ORDNEW.LT.ORDOLD) THEN
          IF (.NOT.TSSEARCH) THEN
            WRITE(LUOUT,1884)
 1884       FORMAT(T3,'@SYMMETRY-F, Descent in symmetry detected.')
            CALL ERREX
          ENDIF
        ELSE
            WRITE(LUOUT,1885)
 1885       FORMAT(T3,'@SYMMETRY-F, Ascent in symmetry detected.')
            CALL ERREX
        ENDIF
        ENDIF
      ENDIF
C
C CHECK TO SEE IF ORIENTATION MATRIX IS EMPTY.  IF IT IS, THEN SET IT EQ
C  TO THE IDENTITY MATRIX.
C
      ZTOTAL=0.D0
      DO 1032 I=1,3
        DO 1033 J=1,3
          ZTOTAL=ZTOTAL+DABS(ORIENT(I,J))
 1033   CONTINUE
 1032 CONTINUE
      IF(ZTOTAL.LT.1.D-8)THEN
        DO 1034 I=1,3
          ORIENT(I,I)=ONE
 1034   CONTINUE
      ENDIF
C
C DEFINE ORIENT TO PUT I AND Ih IN THEIR CANONICAL ORIENTATIONS (C5s
C  ALONG (0,0,1) AND (0,2/SQRT(5),1/SQRT(5)).
C
      IF(FPGRP.EQ.'I h '.OR.FPGRP.EQ.'I   ')THEN
        CALL ZERO(SCRATCH,3*NX)
        SQ5=ONE/DSQRT(5.D0)
        ARGX=DSQRT(0.5D0*(ONE-SQ5))
        ANGL=DACOS(ARGX)
        ANGL2=DACOS(SQ5)
        CALL ROTM(1,ANGL,0,ORIENT)
      ENDIF
C
C DUMP SOME NECESSARIES TO JOBARC
C
      IF(FPGRP(2:2).EQ.'X')ILINEAR=1
#ifdef _DEBUG_LVL0
      Print*, "The Orientation matrices: ORIEN2 and ORIENT"
      Write(6,*)
      CALL OUTPUT(ORIEN2, 1, 3, 1, 3, 3, 3, 0)
      Write(6,*)
      CALL OUTPUT(ORIENT, 1, 3, 1, 3, 3, 3, 0)
      Write(6,*) "The coords in /COORD/ common block"
      Print*, (Q(I), I=1, 3*NATOMS)
      Write(6,*)
      Write(6,*) "The coords NEWQ array "
      Print*, (NEWQ(I), I=1, 3*NATOMS)
      Write(6,*)

#endif
      CALL XGEMM('N','N',3,3,3,ONE,ORIEN2,3,ORIENT,3,ZILCH,
     &     SCRATCH,3)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      CALL OUTPUT(SCRATCH, 1, 3, 1, 3, 3, 3, 0)
#endif

      CALL IPUTREC(20,'JOBARC','LINEAR  ',IONE,ILINEAR)
      RETURN
      END
