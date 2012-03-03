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
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
      SUBROUTINE RETRIEVE(E, V1,V2,STEP,VEC)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "fnamelen.par"
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      LOGICAL XYZIN, NWFINDIF
      CHARACTER*(fnamelen) FNAME
C
#include "cbchar.com"
#include "coord.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $        ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $        XYZTol
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
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
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      CHARACTER*10 CRAP,CRAP2
      DIMENSION V1(NXM6),V2(NXM6*NXM6),STEP(NXM6),VEC(NOPT)
      INTEGER    LuArc
      PARAMETER (LuArc  = 77)
      CHARACTER*(*) ArcFil
      PARAMETER    (ArcFil = 'OPTARC ')
#include "io_units.par"
      DATA CRAP /'***CYCLE  '/
C
      IF( ncycle .NE. 0 )WRITE(6,933)
933   FORMAT(T3,' Retrieving information from last optimization cycle.')
      CALL GFNAME(ARCFIL,FNAME,ILENGTH)
      OPEN(UNIT=LuArc,FILE=fname(1:ilength),STATUS='OLD',
     &     FORM='UNFORMATTED')
C
C NUMBER OF RECORDS IN .OPT IS ONE GREATER THAN THE NUMBER OF
C  OPTIMIZATION CYCLES.
C
      DO 101 K=1,nCYCLE
         READ(LuArc)CRAP2
         IF(CRAP2.NE.CRAP)then
            Write (LuErr,*) ' *Problem with optimization archive file.'
            Call ErrEx
         EndIf
 101  Continue
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "Start Reading Archive file NXM6: ", NXM6
#endif
C
      READ(LuArc)CRAP2,iarchx,ICYCLE,NATOMS,NX,NUNIQUE,NOPT,IPRNT,INR,
     $     IVEC,
     $ IDIE,IMXSTP,ISTCRT,IVIB,ICURVY,ICONTL,IRECAL,INTTYP,IDISFD,
     $ IGRDFD,ICNTYP,ISYM,IBASIS,E,
     $   (V1(I),I=1,NXM6),
     $   (V2(I),I=1,NXM6*NXM6),(STEP(J),J=1,NXM6),(R(J),J=1,NXM6),
     $   (Q(J),J=1,NX),(NCON(I),I=1,NX),(ZSYM(I),I=1,NATOMS),
     &   (VARNAM(I),I=1,NX),(PARNAM(I),I=1,NX),(ISQUASH(J),J=1,NXM6),
     &   (IATNUM(I),I=1,NATOMS),(ATMASS(I),I=1,NATOMS),
     $   (IUNIQUE(J),J=1,NUNIQUE),(NEQ(J),J=1,NXM6),
     $   ((IEQUIV(I,J),J=1,NEQ(IUNIQUE(I))),I=1,NUNIQUE),
     &   (NOPTI(J),J=1,NOPT), (VEC(J),J=1,NOPT),FPGRP,BPGRP,PGRP,
     $   ((ORIENT(I,J),J=1,3),I=1,3)
         IF(CRAP2.NE.CRAP)then
            Write (LuErr,*) ' *Problem with optimization archive file.'
            Call ErrEx
         EndIf
C
#ifdef _DEBUG_LVL0
      Print*, "The COORD COMMON BLOCK AT RETRIEVE BEFORE DECOMPRESS"
      Write(6, "(3F10.5)") (Q (I), I= 1, 3*NATOMS)
      Write(6,*)
      Write(6, "(3F10.5)") (R(I), I= 1, NXM6)
#endif

C
C DECOMPRESS INTERNAL COORINATE VECTOR
C
      XYZIN = IFLAGS(h_IFLAGS_coordinates) .GT. 0
C
      IF (.NOT. XYZIN) THEN 
         DO J=NX,10,-1
            R(J)=R(J-6)
         ENDDO 
         R(8)=R(3)
         R(7)=R(2)
         R(4)=R(1)
         R(1)=0.D0
         R(2)=0.D0
         R(3)=0.D0
         R(5)=0.D0
         R(6)=0.D0
         R(9)=0.D0
       ENDIF
C
#ifdef _DEBUG_LVL0
      Print*, "The COORD COMMON BLOCK AT RETRIEVE"
      Write(6, "(3F10.5)") (Q (I), I= 1, 3*NATOMS)
      Write(6,*)
      Write(6, "(3F10.5)") (R(I), I= 1, NXM6)
#endif

C
      CLOSE(77,STATUS='KEEP')
      RETURN
      END
