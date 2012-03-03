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
      SUBROUTINE ENTRY(BasNam, We_havegeom, Can_do_freq)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      character*(*) BasNam

      CHARACTER*10 TRASH
#include "fnamelen.par"
#include "flags.h"
      CHARACTER*(fnamelen) FNAME
      LOGICAL OLDARC,XYZIN,NWFINDIF,YESNO,IGNORE,OPTRES,
     &        We_havegeom, Can_do_freq
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
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON/RESTART_COM/IGNORE
      COMMON/RESTART2/OPTRES
C 
      INTEGER    LuArc
      PARAMETER (LuArc  = 77)
      CHARACTER*(*) ArcFil
      PARAMETER    (ArcFil = 'OPTARC ')
#include "io_units.par"
C From COMMON:
C     USINT:  IARCH, NCYCLE, NX
C     LUINTS:  Archive unit & name, Abinitio unit & name
C
      Call igetrec(0, 'JOBARC', 'HESSIANM', Length, Scr)
C
C If analytical hessians are available and the Hessian is on 
C the disk we can transfer the logic to simply print out the
C frequencies.
C
      If (Length .GT. 0) Then
         CALL GTFLGS(0,IERR,IPRNT1,INR,ICONTL,IVEC,IDIE,
     &               ICURVY,ISTCRT,IMXSTP,IVIB,  IRECAL,
     &               INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,
     &               IBASIS,IDFGHI,
     &               BasNam)
          If (IVIB .EQ. 1) Then
              Iarch = 1
              Call igetrec(20, 'JOBARC', 'NREALATM', 1, NATOMS)
              Call igetrec(20, 'JOBARC', 'LINEAR  ', 1, ILINEAR)
             If (ILINEAR .EQ. 1) Then
                 NX   = 3*NATOMS
                 NXM6 = NX - 5
             Else
                 NX   = 3*NATOMS
                 NXM6 = NX -6
             Endif
           Endif
           Can_do_freq = (IVIB .GT. 0 .AND. We_havegeom)
           If (Can_do_freq) Then
               Return
           Endif
      Endif

C See if there is an existing archive file
C
      CALL GFNAME(ARCFIL,FNAME,ILENGTH)
      INQUIRE (FILE=FNAME(1:ILENGTH), EXIST=OLDARC)

#ifdef _DEBUG_LVL0
      Print*, "The OPTARC is here?", OLDARC, We_havegeom 
      Print*, "Do we have geom?", We_havegeom
      Print*, "Hessian in JOBARC?", Length
#endif
C
      IF (OLDARC .AND. We_havegeom) Then
C
C This is where prepartion are made to read the optimized geometry. 
C The OPTARC has the final optimized geometry (written by summary.F) 
C and it is read by retrive.F (in geopt) and prepare for the 
C vibrational frequency calculation.
C

         Call igetrec(20, 'JOBARC', 'NREALATM', 1, NATOMS)
         Call igetrec(20, 'JOBARC', 'LINEAR  ', 1, ILINEAR)
         CALL GTFLGS(0,IERR,
     &               IPRNT1,INR,   ICONTL,IVEC,  IDIE,
     &               ICURVY,ISTCRT,IMXSTP,IVIB,  IRECAL,
     &               INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,
     &               IBASIS,IDFGHI,
     &               BasNam)
         IF (ILINEAR .EQ. 1) Then
             NX   = 3*NATOMS
             NXM6 = NX - 5
         ELSE
             NX   = 3*NATOMS
             NXM6 = NX -6 
         ENDIF
         Return
      Endif
C
      IF (OLDARC) THEN
C
C This block is executed during a geometry optimization.
C
         OPEN (LUARC,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &         STATUS='OLD')
         REWIND (LUARC)
C Read until the last recorded cycle

 10      READ (LUARC, END=20) TRASH, IARCH, NCYCLE, NJUNK, NX, NUNIQUE,
     $      NOPT, IPRNT, INR, IVEC, IDIE, IMXSTP, ISTCRT , IVIB ,
     &      ICURVY, ICONTL, IRECAL, INTTYP, IDISFD, IGRDFD, ICNTYP
     &      ISYM, IBASIS
         IF (ICONTL.EQ.0) ICONTL=4
         GOTO 10
 20      CONTINUE
         Close (LuArc, Status='KEEP')
C
C This is a continuation, so set bit indicating existance of
C the archive file. The iarch is tested in geopt to switch the
C the geometry optimization block
C
              idump = 0
              iarch = 1
             itmp   = iflags(18)
         iflags(18) = mod(itmp,100)
         IF (IVIB.NE.0.OR.
     &       IFLAGS(18).EQ. 3.OR.         ! NMR
     &       IFLAGS(18).EQ. 8.OR.         ! J_SO
     &       IFLAGS(18).EQ. 9.OR.         ! J_FC
     &       IFLAGS(18).EQ.10.OR.         ! J_SD
     &       IFLAGS(18).EQ. 4.OR.         ! NMR_SWTCH
     &       IFLAGS(18).EQ. 5    )THEN    ! CONVENTIONAL
            CALL GTFLGS(0,IERR,
     &                  IPRNT1,INR,   ICONTL,IVEC,  IDIE,
     &                  ICURVY,ISTCRT,IMXSTP,IVIB,  IRECAL,
     &                  INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,
     &                  IBASIS,IDFGHI,
     &                  BasNam)
                  itmp = iflags(18)
            iflags(18) = mod(itmp,100)
            XYZIN      = (IFLAGS(68).EQ.1)
         END IF
C
         iflags(18)=itmp
C
C Else for the (IF OLDARC)
C
      ELSE
C
         CALL igetrec(1,'JOBARC','FIRSTRUN',1,iDump)
C
         INQUIRE(FILE=ZFIL,EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(LUOUT,9003)
 9003       FORMAT(T3,'@ENTRY-F, ZMAT file not present.')
            CALL ERREX
         END IF
C
C Endif for The IF (OLDARC)
C 
      END IF
C
C This statement is executed pratically all situations except 
C in the case Freq. followed by an optimizations in which case
c the control is transfered to the first block (see above).
C
      CALL GTFLGS(IDUMP,IERR,IPRNT,INR,ICONTL,IVEC,IDIE,ICURVY,ISTCRT,
     &IMXSTP,IVIB,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     &IDFGHI,BasNam)

      OPTRES = OLDARC.AND.(.NOT.IGNORE)
C
      IF (OPTRES) THEN
         WRITE(LUOUT,9009) NCYCLE
 9009    FORMAT(T3,' JODA restarting optimization with cycle #',I3,'.')
      END IF
C
             itmp = iflags(18)
      iflags(18)  = mod(itmp,100)
      Can_do_freq = (IVIB .GT. 0 .AND. We_havegeom) 
C
      IF((.NOT. Can_do_freq ) .AND.IFLAGS(18).NE.3.
     &   AND.IFLAGS(18).NE.8.AND.IFLAGS(18).NE.9.
     &   AND.IFLAGS(18).NE.10.AND.IFLAGS2(3).NE.1.
     &   AND.IFLAGS(18).NE.4.AND.IFLAGS(18).NE.5.AND.OLDARC)THEN
C
       IF(.NOT.OPTRES) WRITE(LUOUT, 9010) NCYCLE+1
C
      ENDIF
C
      iflags(18) = itmp
C
9010  FORMAT (T3,' JODA beginning optimization cycle #', I3,'.')

      RETURN
      END

