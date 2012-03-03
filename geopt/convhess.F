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
      SUBROUTINE CONVHESS(A,SCRATCH,HC,HI,AT,ID)
C
C     CONVERTS CARTESIAN HESSIAN TO INTERNAL COORDINATE HESSIAN
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "io_units.par"
#include "fnamelen.par"
#include "flags.h"
#include "jodaflags.com"
C
      CHARACTER*(fnamelen) FNAME
      LOGICAL XYZIN, NWFINDIF
      INTEGER TOTNOFBND, TOTNOFANG, TOTNOFDIH
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
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DOUBLE PRECISION A(NX,NXM6),HC(NX,NX),HI(NXM6,NXM6)
      DOUBLE PRECISION AT(NXM6,NX),SCRATCH(NX,NXM6)
C
      IF (ID .GT. 0) THEN
         CALL ZERO(HI, NXM6*NXM6)
         IF (XYZIN) THEN
            IF (iFlags2(h_IFLAGS2_geom_opt) .ge. 3) THEN
               CALL IGETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
               CALL IGETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
               CALL IGETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
               DO IBONDS = 1, TOTNOFBND
                  HI(IBONDS, IBONDS) = 1.0D0
                END DO
                DO IANGS = TOTNOFBND + 1, TOTNOFANG + TOTNOFBND
                  HI(IANGS, IANGS) = 0.25D0
                END DO
                DO IDIHS = TOTNOFBND + TOTNOFANG + 1, TOTNOFANG 
     &                     + TOTNOFBND + TOTNOFDIH
                   HI(IDIHS, IDIHS) = 0.10D0
                END DO
                RETURN
         ELSE IF (iFlags2(h_IFLAGS2_geom_opt) .eq. 2) THEN
                CALL SET_2UNIT_MATRIX(HI, NX)
         ENDIF

        ELSE
            HI(1,1)=1.D0
            HI(2,2)=1.D0
            HI(3,3)=0.25D0
            IF(NXM6.GT.3)THEN
              DO 10 I=4,NXM6-2,3
                 HI(I,I)    =1.D0
                 HI(I+1,I+1)=0.25D0
                 HI(I+2,I+2)=0.1D0
10            CONTINUE
             ENDIF
            RETURN 
        ENDIF
      ENDIF
C
C
C The following block of code is executed only in vibrational
C frequency calculations. The purpose is to bulid a Hessian in
C internal cordinates and write it to the FCMINT file. Note the
C assumption in the transformation: the gradient is zero, hence
C no derivative of B matrix is required. This Hessian can be
C used as a starting Hessian for Transition State searches.
C In the case of frequency calculations with Cartesians let's not do
C any transformations and leave the Hessian in Cartesian
C Coordinates.  Ajith Perera 07/2003.
C
      IF (.NOT.XYZIN) THEN
         CALL MODMATMUL(SCRATCH,HC,A,NX,NX,NXM6,NX,NX,NXM6)
         CALL MTRANSP(A,AT,NX,NXM6,NX,NXM6)
         CALL MODMATMUL(HI,AT,SCRATCH,NXM6,NX,NXM6,NXM6,NX,NXM6)
      END IF
C
C Note that the Cartesian Hessian is already in the JOBARC record
C HESSINAM written by vdint for analytical Hessians and symcor for
C numerical Hessian calculations. Since there can be ordering
C issues (vmol vs ZMAT order) let's write another record. We know
C that this is in ZMAT order. 01/2006, Ajith Perera.
C
      CALL DPUTREC(20,'JOBARC','CART_HES',NX*NX,HC)
C
#ifdef _DEBUG_LVL0
      Print*, "Writing Internal or Cartesian Hessian"
      Print*, "XYZIN, NX, NXM6:", XYZIN, NX, NXM6 
      CALL OUTPUT(HI, 1, NX, 1, NX, NX, NX, 1)
      print *,'Cartesian Hessian ### '
       CALL OUTPUT(HC, 1, NX, 1, NX, NX, NX, 1)
#endif     
      IF(ID.EQ.-1)THEN
         CALL GFNAME('FCMINT  ',FNAME,ILENGTH)
         OPEN(UNIT=50,FILE=FNAME(1:ILENGTH),FORM='FORMATTED',
     &        STATUS='UNKNOWN')
         IF (.NOT.XYZIN) THEN
            WRITE(LUOUT,112)
112         FORMAT(T3,'@CONVHESS: the FCMINT is in internals')
            WRITE(50,'((3(E18.12,1X)))')((HI(I,J),J=1,NXM6),I=1,NXM6)
C
C Write the internal Hessian to JOBARC even though
C we have it in the FCMINT file.
C
            CALL DPUTREC(20,'JOBARC','INTR_HES',NXM6*NXM6,HI)
         ELSE
            WRITE(LUOUT,113)
113         FORMAT(T3,'@CONVHESS: Write Cartesians FCMINT file')
            WRITE(50,'((3(E18.12,1X)))')((HC(I,J),J=1,NX),I=1,NX)
            WRITE(*,*) "NX: ",NX
            WRITE(*,'((3(E18.12,1X)))')((HC(I,J),J=1,NX),I=1,NX)
         ENDIF
          CLOSE(UNIT=50,STATUS='KEEP')
      ENDIF
C
      IF(IPRNT.GE.500)THEN
       WRITE(LUOUT,110)
 110   FORMAT(T3,' @CONVHESS-I, Full Cartesian Hessian: ')
       WRITE(LuOut,90)((I,J,HC(I,J),J=1,I),I=1,NX)
       WRITE(LUOUT,111)
 111   FORMAT(T3,' @CONVHESS-I, Full internal coordinate Hessian: ')
       WRITE(LuOut,90)((I,J,HI(I,J),J=1,I),I=1,NXM6)
 90    FORMAT(3('[',I3,',',I3,']',1X,F9.6,1X),'[',I3,',',I3,']',1X,
     & F9.6)
      ENDIF
C
      RETURN
      END
