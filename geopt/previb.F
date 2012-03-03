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
      SUBROUTINE PREVIB(NRX,H,SCRATCH,GRAD)
C
C FRONT END TO NORMAL COORDINATE CALCULATION.  COMPRESSES DUMMY
C   ATOMS OUT OF VECTORS PASSED THROUGH COMMON (VARIABLE NAMES,
C   ATOMIC MASSES, COORDINATES) AND RETURNS THE NUMBER OF REAL
C   ATOMS * 3.  ALSO CALCULATE VECTOR WHICH RELATES ABSOLUTE
C   Z-MATRIX POSITION WITH "REAL ATOM" POSITION.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
      PARAMETER (TOL = 1.D-3)
      DIMENSION IREPOS(MXATMS)
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
#include "cbchar.com"
#include "coord.com"
C
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
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/OPTCTL/IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     &   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     &   XYZTOL
      COMMON/FLAGS/IFLAGS(100),IFLAGS2(500)
      COMMON/ORIENT/ORIENT(3,3)
#include "io_units.par"
      CHARACTER*5 JUNK(MXATMS)
      DIMENSION SCRATCH(3*NX),H(NX,NX),GRAD(NX)
C
C DETERMINE NUMBER OF REAL ATOMS (NOT DUMMIES).  THIS
C  DETERMINES THE RANK OF THE HESSIAN.  ALSO, SIMULTANEOUSLY
C  REMOVE DUMMY ATOM ENTRIES FROM Q,ATMASS AND ZSYM VECTORS.
C
C
C FIRST GENERATE A COMPRESSED HESSIAN.  REMOVE DUMMY ATOM COLUMNS
C  FIRST. CODED SO THAT THE ARRAY CAN SUBSEQUENTLY BE PASSED INTO
C  VIB1 AS A NRX X NRX ARRAY.
C
      IX=0
      JX=1
      DO 100 J=1,NX
      DO 100 I=1,NX
       INK=1+(I-0.1)/3
       JNK=1+(J-0.1)/3
       IF(ATMASS(INK).GT.0.1.AND.ATMASS(JNK).GT.0.1)THEN
       IX=IX+1
       H(IX,JX)=H(I,J)
       IF(IX.EQ.NX)IX=0
       IF(IX.EQ.0)JX=JX+1
      ENDIF
100   CONTINUE
C
      Call igetrec(20,'JOBARC','NREALATM',1,Nreals)
      NRX = 3*NREALS
C
C The Cartesian Hessian ordered according to the input and without dummy
C atom contributions. Used in Coriolis Coupling modlue in a2proc.
C Ajith Perera, 02/2009
c
      CALL DPUTREC(20,'JOBARC','CARTHESC',NRX*NRX,H)
      CALL ZERO(SCRATCH,3*NX)
C
#ifdef _DEBUG_LVL0
C      Write(6,*) "@-previb, The Hessian"
       CALL output(h, 1, NRX, 1, NRX, NRX, NRX, 1)
       Write(6,*)
       Write(6,*) "@-previb, The griadients"
       Write(6, "(3F10.5)") (GRAD(I), I=1,NX)
#endif
C
C
C NOW COMPRESS THE VARIOUS ARRAYS.
C
 
      NREAL=0
      NRX=0
      CALL ZERO(SCRATCH,3*NX)
      CALL IZERO(IREPOS,NATOMS)
      DO 3 I=1,NATOMS
      IF(ATMASS(I).GE.0.5D0)THEN
       NREAL=NREAL+1
       IREPOS(I)=NREAL
       SCRATCH(NREAL)=ATMASS(I)
       JUNK(NREAL)=ZSYM(I)
       CALL XDCOPY(3,GRAD(3*I-2),1,GRAD(3*NREAL-2),1)
       CALL VADD(SCRATCH(NX+3*NREAL-2),SCRATCH(NX+3*NREAL-2)
     &,Q(3*I-2),3,1.D0)
      ENDIF
 3    CONTINUE
      CALL ZERO(Q,NX)
      CALL ZERO(ATMASS,NATOMS)
      NRX=3*NREAL
      CALL VADD(Q,Q,SCRATCH(NX+1),NRX,1.D0)
      CALL VADD(ATMASS,ATMASS,SCRATCH,NREAL,1.D0)
      DO 10 I=1,NREAL
      ZSYM(I)=JUNK(I)
 10   CONTINUE
      DO 11 I=NREAL+1,NATOMS
      ZSYM(I)='     '
 11   CONTINUE
      IF(NREAL.EQ.0)THEN
       WRITE(LUOUT,8702)
 8702  FORMAT(T3,'@PREVIB-F, No real atoms found in ATMASS string.')
       Call ErrEx
      ENDIF
      CALL IPUTREC(20,'JOBARC','DUMSTRIP',NATOMS,IREPOS)
C
C COMPRESS DUMMY ATOM STUFF OUT OF BMATRIX
C
      IF(IFLAGS(68).NE.1) THEN
C
       CALL DGETREC(20,'JOBARC','BMATRIX ',NX*NXM6,SCRATCH)
       IX=0
       JX=0
       DO 110 I=1,NX
        INK=1+(I-0.1)/3
        JX=JX+1
        IF(IREPOS(INK).NE.0)THEN
         IX=IX+1
         IF(IX.NE.JX)THEN
          IOFF=1+(IX-1)*NXM6
          JOFF=1+(JX-1)*NXM6
          CALL XDCOPY(NXM6,SCRATCH(JOFF),1,SCRATCH(IOFF),1)
         ENDIF
        ENDIF
110    CONTINUE
       CALL DPUTREC(20,'JOBARC','BMATRIXC',3*NREAL*NXM6,SCRATCH)
C
      ENDIF
C
      RETURN
      END
