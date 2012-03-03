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
      SUBROUTINE VIB1(H,HBUF,HINT,SCRATCH,ECKART,BUF1,BMAT,GRAD,
     &                POLDER,POLTMP,RINT,RDPL,NRX,IPROJ,SCR2)
C
C NORMAL MODE ANALYSIS FOR ACES II.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "io_units.par"
#include "fnamelen.par"
C TOLERANCE FOR DEGENERACY OF EIGENVALUES.
      PARAMETER (DEGTLX = 1.D-6)
      INTEGER*8 NIMAGD, NRXD
      CHARACTER*1 XYZ(3),FLAG(3*MXATMS)
      CHARACTER*11 TYPE(3*MXATMS)
      CHARACTER*4 IRREP(3*MXATMS)
      CHARACTER*10 LABEL
      CHARACTER*4 TMPSTR, CIRREP(3*MXATMS), CHA_TMP1,
     &            CHA_TMP2, IRPCHAR(8)
      CHARACTER*(fnamelen) FNAME
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
#include "cbchar.com"
#include "coord.com"
C
      LOGICAL YESNO,DIPOLE,READIS,RAMAN,POLDER_EXIST,
     &        SYMINFO_EXSIST
      INTEGER EOR(8, 8), KA(8, 3)
      DIMENSION SCRATCH(3*NRX),ECKART(2*NRX),IPT(3),BUF1(NXM6,NRX),
     &          SCR2(NRX*NRX),SYOP(72)
      DIMENSION BMAT(NXM6,NRX),GRAD(NRX), POLDER(9*NRX), POLTMP(9,NRX)
      DIMENSION H(NRX,NRX),HBUF(NRX,NRX),HINT(NRX),NORD(6*MXATMS),
     &          RINT(NRX),RDPL(NRX)
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
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
C
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
C
      Character*4 FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)
      DATA XYZ /'X','Y','Z'/
      DATA HALF,FACT01 /0.5D0,2.85914392D-03/
      NREAL=NRX/3
      DEGTOL=DEGTLX
C
      Print*, "The point group", FPGrp, BPGrp, PGrp 
      ione = 1
      ibad = 0
      call igetrec(-1, 'JOBARC', '12SWITCH', ione, ibad)
      if (ibad .eq. 1) then
         write(6,*) ' Atom 3 is to be connected to atom 2'
         write(6,*) ' Internal vibrational coordinates will be wrong'
         write(6,*) ' see Reconstruct Zmat at beginning of output file'
C         if (iflags2(3) .eq. 1) then
C            write(6,*) ' RESRAMAN output is definitely nonsense '
C         call errex
C         endif
      endif
C
C RELAX CRITERION FOR DEGENERACY IF THIS IS A FINITE DIFFERENCE 
C CALCULATION 
C
      RAMAN=(IFLAGS2(151).EQ. 1)
      IF(IVIB.ge.2)DEGTOL=DEGTOL*10.D0
C
C RELAX FURTHER IF IT IS SINGLE-SIDED.
C
      IF(IVIB.GE.2.AND.IGRDFD.EQ.1)DEGTOL=DEGTOL*10.D0
      istop=0
      WRITE(*,*) 'NRX in VIB1:', NRX ! Watson debug
      IF(MOD(NRX,3).NE.0)THEN
       WRITE(6,8701)
 8701 FORMAT(T3,'@VIB1-F, Something very wrong with passed parameters.')
      Call ErrEx
      ENDIF
      IF(IPRNT.GE.8)THEN
       WRITE(LUOUT,8708)NREAL*3
 8708  FORMAT(T3,'@VIB1-I, The rank of the Hessian is ',I3)
      ENDIF
      IF(NREAL.EQ.0)THEN
       WRITE(LUOUT,8702)
 8702  FORMAT(T3,'@VIB1-F, No real atoms found in ATMASS string.')
       Call ErrEx
      ENDIF
C
C READ IN ISOTOPIC MASSES IF REQUIRED
C
      CALL GFNAME('ISOMASS ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=READIS)
      IF(READIS)THEN
       OPEN(UNIT=70,FILE=FNAME(1:ILENGTH),FORM='FORMATTED',STATUS='OLD')
       READ(70,*)(ATMASS(J),J=1,NREAL)
       WRITE(6,3030)
3030   FORMAT(T3,'@VIB1-I, Non-standard isotopic masses used.')
       WRITE(6,3031)
3031   FORMAT(T16,'Atom',T29,'Atomic mass')
       DO 3040 IPOS=1,NREAL
        WRITE(6,3032)IPOS,ATMASS(IPOS)
3032    FORMAT(T17,I2,T30,F9.5)
3040   CONTINUE
       CLOSE(UNIT=70,STATUS='KEEP')
      ENDIF
C
C MASS-WEIGHT HESSIAN MATRIX
      Write(6,*) "Before mass"
      call output(h, 1, nrx, 1, nrx, nrx, nrx, 1)
C
      DO 15 I=1,NRX
      DO 15 J=1,NRX
      FACT=DSQRT(ATMASS(1+(I-1)/3)*ATMASS(1+(J-1)/3))
      IF(FACT.LT.1.D-3)THEN
       H(I,J)=0.D0
      ELSE
       H(I,J)=H(I,J)/FACT
      ENDIF
   15 CONTINUE
      Write(6,*) "The Hess in vib1"
      call output(h, 1, nrx, 1, nrx, nrx, nrx, 1)
C
C WRITE OUT FCM TO SCRATCH FILE IN CASE IT IS NEEDED LATER.
C
C
C FIRST MAKE SURE IT IS NOT THERE.
C
      INQUIRE(FILE='FCMSCR',EXIST=YESNO)
      IF(YESNO)THEN
       INQUIRE(FILE='FCMSCR',OPENED=YESNO,NUMBER=LUOLD)
       IF(YESNO)THEN
        CLOSE(LUOLD,STATUS='DELETE')
       ELSE
        OPEN(UNIT=78,FILE='FCMSCR',STATUS='OLD',FORM='UNFORMATTED')
        CLOSE(UNIT=78,STATUS='DELETE')
       ENDIF
       WRITE(LUOUT,341)
 341   FORMAT(T3,'@VIB1-I, File FCMSCR was on disk and was deleted.')
      ENDIF
      OPEN(UNIT=78,FILE='FCMSCR',STATUS='NEW',FORM=
     &'UNFORMATTED')
C     REWIND(78)
      WRITE(78)((H(I,J),J=1,NRX),I=1,NRX)
C
C DIAGONALIZE - EIGENVALUES AND EIGENVECTORS RETURNED SORTED -
C   EIGENVALUES ARE IN DIAGONAL ELEMENTS OF H, EIGENVECTORS ARE
C   COLUMNS OF HBUF.  EIGENVECTORS ARE IN *MASS-WEIGHTED* CARTESIAN
C   COORDINATES HERE.
C
      CALL EIG(H,HBUF,NRX,NRX,0)
C
C CHECK FOR NEGATIVE EIGENVALUES - SET BIT TO PRINT 'i' IF ENCOUNTERED.
C
      DO 26 I=1,NRX
      FLAG(I)=' '
      IF(H(I,I).LT.0.D0)FLAG(I)='i'
   26 CONTINUE
C
C CONVERT EIGENVALUES TO CORRESPONDING FREQUENCIES
C
      DO 25 I=1,NRX
   25 H(I,I)=DSQRT(DABS(H(I,I)))*5.14048D03
C
C WRITE NORMCO FILE
C (used to produce interfaces to various visualization programs such as
C  MOLDEN or HYPERCHEM.)
C
      OPEN(UNIT=56,FILE='NORMCO',FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(56,7500)
7500  FORMAT('% mass weighted coordinates')
      IOFF=1
      DO 20 I=1,NRX/3
       WRITE(56,'(3F20.10)')(SQRT(ATMASS(I))*Q(J),J=IOFF,IOFF+2)
       IOFF=IOFF+3
20    CONTINUE
      DO 21 IMODE=1,NRX
       WRITE(56,7501)
7501   FORMAT('% frequency')
       WRITE(56,'(F20.10,A1)')H(IMODE,IMODE),FLAG(IMODE)
       WRITE(56,7502)
7502   FORMAT('% normal coordinate ')
       WRITE(56,'((3F20.10))')(HBUF(J,IMODE),J=1,NRX)
21    CONTINUE
C Note that the intensities are written to the NORMCO file below.
      CALL XDCOPY(NRX,H,NRX+1,SCRATCH,1)
      CALL DPUTREC(20,'JOBARC','FREQUENC',NRX,SCRATCH)
C
C CHECK EIGENVECTORS AGAINST ECKART CONDITIONS.  ONLY DO THIS IF
C FREQUENCIES ARE UNPROJECTED.
C
      IF(IVIB.NE.3.OR.(IVIB.EQ.3.AND.IFLAGS(80).NE.0))THEN
C
CJDW 9/11/96. Note that CHKECK makes dimension of HBUF MASS * DISTANCE.
C             Below (in this part of the IF-block) we divide by the square
C             root of mass to recover the original eigenvectors. It is
C             these MASS**1/2 * DISTANCE eigenvectors which are printed
C             out under "Normal Coordinates".
C
       CALL CHKECK(NRX,HBUF,ECKART)
C
C NOW SORT ECKART VECTOR SO THAT WE CAN FIND THE MOST "ROTATIONLIKE" AND
C  "TRANSLATIONLIKE" MODES.  NEED TO PUT IN A CHECK HERE FOR LINEARITY
C  OF THE MOLECULE.
C
C INITIALIZE ALL TYPES TO VIBRATION
C
       DO 1032 I=1,NRX
        TYPE(I)=' VIBRATION '
 1032  CONTINUE
C
C PREPARE VECTOR WHICH WILL EVENTUALLY CONTAIN POINTERS TO THE TRANS.
C   AND ROTS.
C

       DO 5702 I=1,NRX
        NORD(I)=I
        NORD(NRX+I)=I
 5702  CONTINUE
C
C TRANSLATIONS FIRST.  ALWAYS THREE OF THESE.
C
       CALL PIKSR2(NRX,ECKART(NRX+1),NORD(NRX+1))
       DO 5703 I=0,2
       TYPE(NORD(2*NRX-I))='TRANSLATION'
 5703  CONTINUE
C
C CHECK INERTIA TENSOR FOR ZERO ELEMENT, MEANING MOLECULE IS LINEAR.
C
       CALL INERTIA(SCRATCH,.FALSE.)
       ILINER=0
       DO 5704 I=1,9,4
       IF(DABS(SCRATCH(I)).LT.1.D-5)ILINER=ILINER+1
 5704  CONTINUE
       ZIX=SCRATCH(1)
       ZIY=SCRATCH(5)
       ZIZ=SCRATCH(9)
       IF(ILINER.GT.1)THEN
        WRITE(LUOUT,5705)
 5705   FORMAT(T3,'@VIB1-F, Problem with inertia tensor.')
        WRITE(LUOUT,*)(SCRATCH(I),I=1,9)
        Call ErrEx
       ENDIF
C
C NOW ASSIGN ROTATIONAL MODES.
C
       CALL PIKSR2(NRX,ECKART,NORD)
       DO 5706 I=0,2-ILINER
       TYPE(NORD(NRX-I))=' ROTATION '
 5706  CONTINUE
C
C PUT EIGENVECTORS BACK INTO SQRT(MASS)*DISTANCE DIMENSION
C
       DO 323 I=1,NRX
       DO 333 J=1,NRX
       HBuf(J,I)=HBuf(J,I)/DSQRT(ATMASS(1+(J-1)/3))
 333   CONTINUE
       CALL NORMAL(HBUF(1,I),NRX)
 323   CONTINUE
      ELSE
       DO 334 I=1,NRX
        TYPE(I)=' VIBRATION '
        IF(H(I,I).LT.0.01)TYPE(I)=' --------- '
334    CONTINUE
      ENDIF
C
C NOW CHECK TO SEE IF THE DIPOLE DERIVATIVES ARE AVAILABLE.  IF
C  THEY ARE, THEN TRANSFORM THEM TO NORMAL COORDINATES AND DETERMINE
C  THE INTENSITY OF THE BANDS.
C
      INQUIRE(FILE='DIPDER',EXIST=DIPOLE)
      IF(DIPOLE)THEN
C
C DIPOLE DERIVATIVES.
C
C ANALYTIC CALCULATION
C
       IF(IVIB.LE.1)THEN
C
C GET INFORMATION ABOUT ATOMIC ORDERING
C
        CALL IGETREC(20,'JOBARC','DUMSTRIP',NATOMS,NORD)
        CALL IGETREC(20,'JOBARC','MAP2ZMAT',NATOMS,NORD(NATOMS+1))
        OPEN(UNIT=90,FILE='DIPDER',FORM='FORMATTED',STATUS='OLD')
C
C READ THE DIPDERS INTO THE SCRATCH VECTOR AND MASS WEIGHT THEM.
C
        JBOT=1
        DO 1033 I=1,3
         READ(90,*)
         JBOT0=1+(I-1)*3*NREAL
         DO 1034 J=1,NREAL
          JZ=NORD(NATOMS+J)
          JZREAL=NORD(JZ)
          JBOT=JBOT0+(JZREAL-1)*3
          READ(90,'(4F20.10)')ZJUNK,(SCRATCH(K),K=JBOT,JBOT+2)
          Z=1.D0/DSQRT(ATMASS(JZREAL))
         DO 1049 K=JBOT,JBOT+2
           SCRATCH(K)=SCRATCH(K)*Z
1049      CONTINUE
1034     CONTINUE
1033    CONTINUE
        CLOSE(UNIT=90,STATUS='KEEP')
       ELSE
          OPEN(UNIT=90,FILE='DIPDER',FORM='FORMATTED',STATUS='OLD')
          JBOT=1
          DO 1043 I=1,3
             READ(90,*)
             DO 1044 J=1,NREAL
                READ(90,'(4F20.10)')ZJUNK,(SCRATCH(K),K=JBOT,JBOT+2)
                Z=1.D0/DSQRT(ATMASS(J))
                DO 1059 K=JBOT,JBOT+2
                   SCRATCH(K)=SCRATCH(K)*Z
 1059           CONTINUE
                JBOT=JBOT+3
 1044        CONTINUE
 1043     CONTINUE
C
C Read The Polarizability Derivative file, Ajith and John 09/98.
C
        INQUIRE(FILE='POLDER',EXIST=POLDER_EXIST)
        IF (POLDER_EXIST) THEN
           OPEN(UNIT=91,FILE='POLDER',FORM='FORMATTED',STATUS='OLD')
           JBOT=1
           DO 1045 I=1,9
              READ(91,*)
              DO 1046 J=1,NREAL
                 READ(91,'(4F20.10)')ZJUNK,(POLDER(K),K=JBOT,JBOT+2)
                 Z=1.D0/DSQRT(ATMASS(J))
                 DO 1047 K=JBOT,JBOT+2
                    POLDER(K)=POLDER(K)*Z
 1047            CONTINUE
                 JBOT=JBOT+3
 1046         CONTINUE
 1045      CONTINUE
C
           CLOSE(UNIT=91,STATUS='KEEP')
        ENDIF
C
        CLOSE(UNIT=90,STATUS='KEEP')
       ENDIF
C
C TRANSFORM THEM TO THE NORMAL COORDINATE REP. AND WRITE TO SCRATCH FILE
C
C       WRITE(6,*) "POLDER just after Read"
C       CALL OUTPUT(POLDER, 1, 9, 1, NRX, 9, 1)
       CALL MTRAN2(HBUF,NRX)
       CALL MODMATMUL(SCRATCH,HBUF,SCRATCH,NRX,NRX,3,NRX,NRX,3)
       CALL MODMATMUL(POLDER,HBUF,POLDER,NRX,NRX,9,NRX,NRX,9)
       CALL MTRAN2(HBUF,NRX)
C       CALL OUTPUT(SCRATCH, 1, 3, 1, NRX, 3, NRX, 1)
C       CALL OUTPUT(POLDER, 1, 9, 1, NRX, 9, NRX, 1)
C
C Reorder The Polarizability Derivative Matrix and Diagonalize
C Ajith and John 09/98
C
       IF (POLDER_EXIST) THEN
          DO I= 1, NRX
             DO J= 1, 9
                POLTMP(J, I) = POLDER(NRX*(J-1) + I)
             END DO
          END DO
C          CALL OUTPUT(POLDER, 1, 9, 1, NRX, 9, NRX, 1)
          DO I= 1, NRX
             CALL EIG(POLTMP(1,I), POLDER, 3, 3, 1)
          END DO
C          CALL OUTPUT(POLTMP, 1, 9, 1, NRX, 9, NRX, 1)
       END IF
       OPEN(UNIT=95,FILE='DIPSCR',FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND(95)
       WRITE(95)SCRATCH
       CLOSE(95)
C
C Calculate Raman Intensities and Depolarization Ratios.
C Ajith and John, 09/98
C
       IF (POLDER_EXIST) THEN
          DO IMODE = 1, NRX
             ALPHA = (1.0D0/3.0D0)*(   POLTMP(1,IMODE)
     &                               + POLTMP(5,IMODE)
     &                               + POLTMP(9,IMODE))
             BETAS = (1.0D0/2.0D0)*(
     &               (POLTMP(1,IMODE)-POLTMP(5,IMODE))**2
     &             + (POLTMP(5,IMODE)-POLTMP(9,IMODE))**2
     &             + (POLTMP(9,IMODE)-POLTMP(1,IMODE))**2)
             RINT(IMODE) = 0.0784159D0*(45.0D0*ALPHA**2 + 7.0D0*BETAS)
             IF (ALPHA .NE. 0.0D0 .AND. BETAS .NE. 0.0D0) THEN
                RDPL(IMODE) = 3.D0*BETAS/(45.0D0*ALPHA**2 + 4.0D0*BETAS)
             ELSE IF (BETAS .EQ. 0.0D0) THEN
                RDPL(IMODE) = 0.00D0
             ELSE
                RDPL(IMODE) = 0.75D0
             ENDIF
             IF (H(IMODE,IMODE) .LT. 1.0D-3) RDPL(IMODE) = 0.0D0
          END DO
       END IF
C
C COMPUTE INDIVIDUAL BAND INTENSITIES IN KM/MOL.
C
       CALL ZERO(HINT,NRX)
       DO 1035 I=1,NRX
        DO 1036 J=1,3
         Z=SCRATCH(NRX*(J-1)+I)
         HINT(I)=HINT(I)+Z*Z
1036    CONTINUE
1035   CONTINUE
       DO 1038 I=1,NRX
        HINT(I)=HINT(I)*974.868
c     o write the intensities to the NORMCO file
        WRITE(56,7503)
7503    FORMAT('% IR intensity')
        WRITE(56,'(F20.10)')HINT(I)
1038   CONTINUE
      ENDIF
C
C COMMENCE PRINTING OUT THE STUFF, FREQUENCIES FIRST AND THEN
C  NORMAL COORDINATES.
C
 900   IBOT=1
C
C PUT MOLECULE AND EIGENVECTORS IN PA ORIENTATION, SO THAT THE IRREPS
C  CAN BE DETERMINED.
C
C
      CALL STDORT(Q,SCRATCH,3*NREAL,0)

#ifdef _DEBUG_LVLM1
       Write(6,*)
       Write(6,*) "@-VIB1 The normal modes in stand. orientation"
       call output (HBUF, 1, NRX, 1, NRX, NRX, NRX, 1)
#endif
      DO 5002 I=IBOT,NRX
       CALL STDORT(HBUF(1,I),SCRATCH,3*NREAL,0)
 5002 CONTINUE
#ifdef _DEBUG_LVLM1
       Write(6,*)
       Write(6,*) "@-VIB1 The normal modes in sym. orientation"
       call output (HBUF, 1, NRX, 1, NRX, NRX, NRX, 1)
#endif
C
C  DEGENERATATE "Computational IRREPS" FOR FREQUENCIES NOW.
C
c
       CALL ZERO(SCRATCH,3*NRX)
C
C START AT BOTTOM OF LIST.
C
       I=IBOT
C
C CHECK HERE FOR LEVEL OF DEGENERACY AND PASS ALL MEMBERS OF A
C  DEGENERATE SET TO THE GENERATOR
C
        DO 6215 J=IBOT,NRX
        IRREP(J)='----'
 6215   CONTINUE
C
C
C The follwing block of code is completely replaced. The reason
C is that some one other than I (or under my watch) have decided
C to limit the symmetry assignment to the computational point group
C rather than the full point. This is fixed and the logic is simplified.
C For the time being I am keeping this block of code but at a latter
C stage when VIB1.F is cleaned up, all of this will be removed. 
C Ajith Perera, 006/2007.

CSSS          if(ndeg.gt.1)then
C
C RESOLVE DEGENERATE NORMAL COORDINATES ACCORDING TO PURE IRREDUCIBLE
C REPRESENTATIONS OF ABELIAN SUBGROUP
C
CSSS           call igetrec(20,'JOBARC','COMPORDR',1,iorcomp)
CSSS           call dputrec(20,'JOBARC','HOLDDEGQ',5*nrx,hbuf(1,i))
CSSS           call resolve(ndeg,hbuf(1,i),scratch,syop,nord,
CSSS     &                  eckart,scr2,nrx,nrx/3,natoms,iorcomp,ierr)
CSSS          endif
CSSS
cSSS           if(ierr.eq.1)then
cSSS            write(6,*)' normal coordinate cannot be resolved '
cSSS            call dgetrec(20,'JOBARC','HOLDDEGQ',ndeg*nrx,
cSSS     &                  hbuf(1,i))
cSSS           endif
cSSS
C I have done some clean up here. IRREPS() needs to be called for both
C degenerate and non-degenerate cases with the correct arguments to Hbuf
C and CIrrep (controlled by J).
C 01/2006, Ajith Perera.
CSSS          do j = 0, ndeg - 1
CSSS          CALL IRREPS(Q,PGrp,1,HBuf(1,I+J),NRX,CIrrep(i+J),Scratch)
CSSS          enddo
CSSS           do j = 0, ndeg-1
CSSS           CALL IRREPS(Q,FPGrp,NDEG,HBuf(1,I),NRX,Irrep(i),Scratch)
CSSS            do j= 0, ndeg-1
CSSS              Print*, TMPSTR, Irrep(i+J)
CSSS              TMPSTR = IRREP(I)
CSSS              IRREP(I+J) = TMPSTR
CSSS              CIRREP(I+J) = TMPSTR
CSSS           end do
CSSS           Print*, (IRREP(J), J=1, 10)
CSSS          else 
CSSS             CALL IRREPS(Q,PGrp,1,HBuf(1,I),NRX,Irrep(i),Scratch) 
CSSS             CIRREP(I) = IRREP(I) 
CSSS          endif 
CSSS
C ENDIF for the (TYPE(I).EQ.' VIBRATION ') block.
CSSS        ENDIF
CSSS
CSSS        I=I+NDEG
CSSS       IF(I.LE.NRX)GOTO 6218
C
 6218  NDEG = 1
C
      IF (TYPE(I).EQ.' VIBRATION ')THEN
         DO J=1,4
            P=1.D0-H(I+J,I+J)/H(I,I)
            IF(DABS(P).LT.DEGTOL)NDEG=NDEG+1
         ENDDO
 
         if (ndeg.gt.1)then
C 
C RESOLVE DEGENERATE NORMAL COORDINATES ACCORDING TO PURE IRREDUCIBLE
C REPRESENTATIONS OF ABELIAN SUBGROUP
C
           call igetrec(20,'JOBARC','COMPORDR',1,iorcomp)
           call dputrec(20,'JOBARC','HOLDDEGQ',5*nrx,hbuf(1,i))
           call resolve(ndeg,hbuf(1,i),scratch,syop,nord,
     &                  eckart,scr2,nrx,nrx/3,natoms,iorcomp,ierr)
         endif       
C
         CALL IRREPS(Q,FPGrp,NDEG,HBuf(1,I),NRX,Irrep(i),Scratch)
         DO J = 1, NDEG - 1
                  TMPSTR = IRREP(I)
            IRREP(I + J) = TMPSTR
         ENDDO
      ENDIF 
      I = I + NDEG
      IF (I .LE. NRX) GO TO 6218 

C
C Isn't what's below a repeat from above? Someone has done something
C crazy. 01/2006, Ajith Perera.
C
      IF (.FALSE.) THEN
C
C DEGENERATATE "Computational IRREPS" FOR FREQUENCIES NOW.
C
       CALL ZERO(SCRATCH,3*NRX)

C START AT BOTTOM OF LIST.
C
       I=IBOT
C
        DO 2215 J=IBOT,NRX
        IRREP(J)='----'
 2215   CONTINUE
 2218    NDEG=1
        IF(TYPE(I).EQ.' VIBRATION ')THEN
         DO 3215 J=1,4
          P=1.D0-H(I+J,I+J)/H(I,I)
          IF(DABS(P).LT.DEGTOL)NDEG=NDEG+1
 3215 CONTINUE
          IF(TYPE(I).EQ.' VIBRATION ')
     &    CALL IRREPS(Q,FPGrp,NDeg,HBuf(1,I),NRX,Irrep(i),Scratch)
         DO 3216 J=1,NDEG-1
C These two lines were not in the block above. Let's copy them there.
          TMPSTR=IRREP(I)
          IRREP(I+J)=TMPSTR
 3216    CONTINUE
        ENDIF
        I=I+NDEG
        IF(I.LE.NRX)GOTO 2218
      END IF
C
C Appropriate averaging for degenerate irreps, Ajith and John 09/98
C
      IOFF = 1
      DO 1399 IMODE = 1, NRX
C
         DO 1499 I = 1, 4
            IF(IRREP(IOFF)(I:I).EQ.'A'.OR.IRREP(IOFF)(I:I).EQ.'B') THEN
               IOFF = IOFF + 1
            ELSE IF (IRREP(IOFF)(I:I) .EQ. 'E') THEN
               RINT(IOFF) = (RINT(IOFF)+RINT(IOFF+1))/4.0D0
               RINT(IOFF+1) = RINT(IOFF)
               IOFF = IOFF + 2
            ELSE IF (IRREP(IOFF)(I:I) .EQ. 'T') THEN
               RINT(IOFF)=(RINT(IOFF)+RINT(IOFF+1)+RINT(IOFF+2))/9.0D0
               RINT(IOFF+1) = RINT(IOFF)
               RINT(IOFF+2) = RINT(IOFF)
               IOFF = IOFF + 3
            ENDIF
 1499    CONTINUE
         IOFF = IOFF + 1
C
 1399 CONTINUE
C
C NOW PUT EVERYTHING BACK INTO ACES2 ORIENTATION.
C
#ifdef _DEBUG_LVLM1
       Write(6,*)
       Write(6,*) "@-VIB1 The normal modes back in stand. orientation"
       call output (HBUF, 1, NRX, 1, NRX, NRX, NRX, 1)
#endif
      CALL STDORT(Q,SCRATCH,3*NREAL,1)
      DO 5001 I=IBOT,NRX
       CALL STDORT(HBUF(1,I),SCRATCH,3*NREAL,1)
 5001 CONTINUE
C
C This record is specifically added for Coriolis coupling constant
C work. Unlike what is in NORMCO file these are symmetry adapted
C and matches with what is on the output.  Ajith Perera, 02/2009.
C
       CALL DPUTREC(20, "JOBARC", "SYM_NMDS", NRX*NRX, HBUF)
       CALL PUTCREC(20, "JOBARC", "VIB_SYMS", NRX*4, IRREP)
C
#ifdef _DEBUG_LVLM1
       Write(6,*)
       Write(6,*) "@-VIB1 The normal modes back in stand. orientation"
       call output (HBUF, 1, NRX, 1, NRX, NRX, NRX, 1)
#endif
C
      WRITE(LUOUT,1000)
1000  FORMAT(T22,'Normal Coordinate Analysis')
C
      IF (RAMAN) THEN 
          WRITE(LUOUT,2222)
          WRITE(LUOUT,2001)
          WRITE(LUOUT,2002)
          WRITE(LUOUT,2222)
          WRITE(LUOUT,2003)
          WRITE(LUOUT,2222)
C
2001  FORMAT(T2,' Irreducible',T18,' Harmonic ',T34,' Infrared ',
     &T48,'Raman', T63, 'Depolarization',T80,' Type ')
2002  FORMAT(T2,'Representation',T18,' Frequency ',T34,
     &' Intensity ',T48, 'Intensity', T63, 'Ratio')
2003  FORMAT(T20,'(cm-1)',T36,'(km/mol)', T48, '(A**4/AMU)')
C
      ELSE 
C
          WRITE(LUOUT,1111)
          WRITE(LUOUT,3001)
          WRITE(LUOUT,3002)
          WRITE(LUOUT,1111)
          WRITE(LUOUT,3003)
          WRITE(LUOUT,1111)
C
3001  FORMAT(T4,' Irreducible',T24,' Harmonic ',T39,' Infrared ',
     &T55,' Type ')
3002  FORMAT(T4,'Representation',T24,' Frequency ',T39,
     &' Intensity ')
3003  FORMAT(T26,'(cm-1)',T40,'(km/mol)')
C   
      ENDIF 
C
      DO 45 I=IBOT,NRX
         IF (RAMAN) THEN
             WRITE(LUOUT,4000)IRREP(I),H(I,I),FLAG(I),
     &                        HINT(I),RINT(I), RDPL(I), 
     &                        TYPE(I)
4000  FORMAT(T7,A4,T15,F10.4,A1,T30,F10.4,T45,F10.4,T60,F10.4,T75,A11)
C
         ELSE
             WRITE(LUOUT,4005)IRREP(I),H(I,I),FLAG(I),HINT(I),TYPE(I)
 4005 FORMAT(T7,A4,T24,F10.4,A1,T39,F10.4,T55,A11)
C
         ENDIF
C
45    CONTINUE
C
      IF (RAMAN) THEN
          WRITE(LUOUT,2222)
      ELSE
          WRITE(LUOUT,1111)
      ENDIF
C
1111   FORMAT(64('-'))
2222   FORMAT(80('-'))
C
      WRITE(LUOUT,*)
C
      WRITE(LUOUT,5000)
5000  FORMAT(T35,' Normal Coordinates ',/,
     &       T26,' [Dimensions are Mass**-1/2 Distance] ')
C
C NOW SET UP POINTER LIST SO THAT TRANSLATIONS AND ROTATIONS ARE
C   NOT PRINTED
C
C8000  ITOP=MIN(IBOT+2,NRX)
 8000 J=IBOT
      IX=1
 4007 IF(TYPE(J).EQ.' VIBRATION ')THEN
       IPT(IX)=J
       IX=IX+1
       J=J+1
       IF(IX.LT.4.AND.J.LT.NRX+1)GOTO 4007
      ELSE
       J=J+1
       IF(J.LT.NRX+1)GOTO 4007
      ENDIF
      IBOT=J
      WRITE(LUOUT,1100)(IRREP(IPT(I)),I=1,IX-1)
      WRITE(LUOUT,1101)(H(IPT(I),IPT(I)),FLAG(IPT(I)),I=1,IX-1)
      WRITE(LUOUT,1102)(TYPE(IPT(I)),I=1,IX-1)
1100  FORMAT(/,14X,3(A4,22X))
1101  FORMAT(14X,3(F7.2,A1,18X))
1102  FORMAT(12X,3(A11,15X))
      WRITE(LUOUT,1103)((XYZ(J),J=1,3),I=1,IX-1)
1103  FORMAT(T10,3(A1,7X),T36,3(A1,7X),T62,3(A1,7X))
       DO 65 JBT=1,NRX-2,3
C
C SKIP OVER UNINTERESTING DUMMY ATOMS.
C
       WRITE(LUOUT,1105)ZSYM(1+(JBT-1)/3),
     &((HBUF(J,IPT(I)),J=JBT,JBT+2),I=1,IX-1)
1105   FORMAT(T1,A,   F7.4,T14,F7.4,T22,F7.4,T32,F7.4,T40,F7.4,
     &T48,F7.4,T58,F7.4,T66,F7.4,T74,F7.4)
c1105   FORMAT(T2,A,1X,F7.4,T14,F7.4,T22,F7.4,T32,F7.4,T40,F7.4,
c     &T48,F7.4,T58,F7.4,T66,F7.4,T74,F7.4)
65     CONTINUE
C       IBOT=IBOT+3
c      IF(IBOT.GE.NRX)GOTO 9000
CJDW 9/11/96. I think we want the following for linear molecules to
C             work. Format statement 1105 also changed.
       IF(IBOT.GT.NRX)GOTO 9000
       GOTO 8000
9000   CONTINUE
C
C CONVERT GRADIENT VECTOR TO NORMAL COORDINATE REPRESENTATION UNLESS
C VIB = FINDIF
C
       IF(IVIB.LT.2.AND.IFLAGS2(103).NE.0)THEN
        WRITE(6,8004)
        WRITE(6,8001)
        WRITE(6,8002)
        WRITE(6,8001)
        JBOT=1
        DO 2044 J=1,NREAL
         Z=1.D0/DSQRT(ATMASS(J))
         DO 2059 K=JBOT,JBOT+2
          GRAD(K)=GRAD(K)*Z
2059     CONTINUE
         JBOT=JBOT+3
2044    CONTINUE
c     
        CALL XGEMM('N','N',1,NRX,NRX,1.0D0,GRAD,1,HBUF,NRX,
     &             0.0D0,SCRATCH,1)
c
c
c Comments on unit conversion  (from JFS)
c   V = 1/2 k x^2 = 1/2 [4 * pi * pi * nu * nu] * mu * x^2
c
c    = 1/2 [ 4 * pi *pi *nu *nu ] Q^2
c
c   or
c
c     V/hc  =  1/2 omega q^2
c
c     where q = sqrt(2 * pi * c * omega /hbar) * Q
c    and omega is the vibrational frequency in cm-1
c
c  therefore, the energy gradient is
c
c    dE/dq = dE/dQ * dQ/dq
c
c      = dE/dQ * sqrt(hbar/2*pi*c*omega)
c
c   ->   FRED = 0.09114 amu^1/2 bohr divided by sqrt(2.997925d10)
c  (the factor c is taken in the sqrt, see below)
c
c      = 5.265D-7
c
        JBOT=1
        DO 3044 J=1,NREAL
         Z=DSQRT(ATMASS(J))
         DO 3059 K=JBOT,JBOT+2
          GRAD(K)=GRAD(K)*Z
3059     CONTINUE
         JBOT=JBOT+3
3044    CONTINUE
C
        SMAX=0.0D0
        DO 8499 I=1,NRX
         IF(TYPE(I).EQ.' VIBRATION ')THEN
          SCRATCH(NRX+I)=(SCRATCH(I)*SCRATCH(I))/H(I,I)
         ENDIF
8499    CONTINUE
        ILOC=ISAMAX(NRX,SCRATCH(NRX+1),1)
        CALL DSCAL(NRX,1.0D0/SCRATCH(NRX+ILOC),SCRATCH(NRX+1),1)
        DO 8500 I=1,NRX
         IF(TYPE(I).EQ.' VIBRATION ')THEN
          FRED=5.26353359D-07
          C=2.997925D10
          GRADRED1=219474.1D0*SCRATCH(I)/(FRED*DSQRT(C*DABS(H(I,I))))
          GRADRED2=27.2116D0*SCRATCH(I)/(FRED*DSQRT(C*DABS(H(I,I))))
          WRITE(6,8003)I,H(I,I),FLAG(I),SCRATCH(I),GRADRED1,
     &                 GRADRED2,SCRATCH(NRX+I)
         ENDIF
8500    CONTINUE
C 
        WRITE(6,8001)
C
8001   FORMAT(T10,58('-'))
8002   FORMAT(T4,'i',T9,'W(I)',T18,'dE/dQ(i)',T31,'dE/dq',T42,'dE/dq',
     &        T52,'|dE/dq|^2',/,
     &        T31,'(cm-1)',T43,'(eV)',T55,'(relative)')
8003   FORMAT(T3,I3,T7,F8.2,A1,T15,F14.10,T30,F10.3,T41,F10.6,
     &        T52,F14.10)
8004   FORMAT(T13,'Gradient vector in normal coordinate representation')
C
cmn Print out second table in which units / symmetries are displayed
cmn for use with xsim (simulation of vibronic coupling issues.
c
       ione = 1
       ifinal = 0
       call igetrec(-1,'JOBARC', 'FINALSYM', ione, ifinal)
c
c to determine symmetries read in irpchar from file syminfo
c
        inquire(file='SYMINFO', exist=SYMINFO_EXSIST)
        if (SYMINFO_EXSIST) then
        open(unit=10, FILE='SYMINFO',STATUS='OLD')
        do i = 1, 8
           read(10,579) j,  (irpchar(i)(k:k),k=1,4)
        enddo
 579    format(i4,2x,4a1)
        close(unit=10,STATUS='KEEP')
        endif
C
c Construct group multiplication table as in vmol
c
C     Determine KA(I,J)
c
      DO I=1,8
        IS=0
        DO K=1,3
          KA(I,K)=(I-1)/(2**(K-1))-2*((I-1)/(2**K))
          IS=IS+KA(I,K)
       enddo
       enddo
c
      DO I=1,8
        DO J=1,8
        IAN=1
        IOR=1
        DO 102 Kk=1,3
          IAN=IAN+MIN0(KA(I,Kk),KA(J,Kk))*2**(Kk-1)
  102   IOR=IOR+MAX0(KA(I,Kk),KA(J,Kk))*2**(Kk-1)
        EOR(I,J)=IOR-IAN+1
      enddo
      enddo
c
      ione = 1
        call igetrec(-1, 'JOBARC', 'COMPNIRR', ione, nirrep)
cYAU - We don't have to keep printing this matrix out.
c      WRITE (6,2112)
c2112  format(/,'   Group multiplication table :',/)
c      DO 1113 I=1,nirrep
c1113  WRITE (6,1114) (EOR(I,J),J=1,nirrep)
c1114  FORMAT (6x,8I4)
c      write(6,*)

 7001  FORMAT(T3,75('-'))
 7004  format(t13, 'Gradient vector data in eV plus symmetry info')
 7007  format(t13, 'Computational Irrep of electronic state ', a4,
     $      ' -> [',i1,']')
 7005  format(t4, 'Irreps of vibrations ',t32,'Frequencies',
     $ t51, 'Gradient', t62, 'Coupling Sym')
 7002  format(t5,'#', t8, 'FULL  COMP  [#]', t30, 'cm-1',t42,'eV', t51,
     $      'dE/dq(i)',t63,'Irrep',t70,'[#]')
 7006  format(t3,i3,t8,A4,t14,A4, '  [',i1,']', t25, F10.2,A1,' ',
     $ F10.7,A1,t48,F10.7,t62,A4, t70,'[',i1,']')

        WRITE(6,7004)
        if (ifinal .eq. 0) then
           CHA_TMP2 = '----'
        else
           CHA_TMP2 = irpchar(ifinal)
        endif
                                                       
       WRITE(6,7007) CHA_TMP2, ifinal
        write(6,*)
        WRITE(6,7001)
        WRITE(6,7005)
        WRITE(6,7002)
        WRITE(6,7001)
c
        JBOT=1
        DO J=1,NREAL
           Z=1.D0/DSQRT(ATMASS(J))
           DO K=JBOT,JBOT+2
              GRAD(K)=GRAD(K)*Z
           ENDDO
           JBOT=JBOT+3
        ENDDO
c
        CALL XGEMM('N','N',1,NRX,NRX,1.0D0,GRAD,1,HBUF,NRX,
     &             0.0D0,SCRATCH,1)
c
        JBOT=1
        DO J=1,NREAL
         Z=DSQRT(ATMASS(J))
         DO K=JBOT,JBOT+2
            GRAD(K)=GRAD(K)*Z
         ENDDO
            JBOT=JBOT+3
         ENDDO
c
        SMAX=0.0D0
        DO I=1,NRX
           IF(TYPE(I).EQ.' VIBRATION ')THEN
              SCRATCH(NRX+I)=(SCRATCH(I)*SCRATCH(I))/H(I,I)
           ENDIF
        ENDDO
C
        DO I=1,NRX
         IF(TYPE(I).EQ.' VIBRATION ')THEN
c
          FRED=5.26353359D-07
          C=2.997925D10
          GRADRED2=27.2116D0*SCRATCH(I)/(FRED*DSQRT(C*DABS(H(I,I))))
          VIBEV=H(I,I) * 27.2116D0 / 219474.1D0
          i1=0
          do k=1, 8
             if (cirrep(i) .eq. irpchar(k)) i1=k
          enddo
          if (i1 .eq. 0) then
             i2 = 0
          else
          I2=EOR(i1, ifinal)
          endif
          if (i2 .ne. 0) then
             CHA_TMP1 = irpchar(i2)
           else
             CHA_TMP1 = '----'
          endif
          WRITE(6,7006)I,IRREP(I),CIRREP(I),I1 , H(I,I),
     $         flag(i), VIBEV,flag(i), GRADRED2, CHA_TMP1,
     &         I2
         ENDIF
         ENDDO
        WRITE(6,7001)

        ENDIF  
C
C CONVERT VIBRATIONAL COORDINATES TO INTERNALS
C
C
C FIRST PUT EIGENVECTORS BACK INTO DISTANCE DIMENSION
C
       DO 623 I=1,NRX
        IF(TYPE(I).EQ.' VIBRATION ')THEN
         DO 633 J=1,NRX
          HBuf(J,I)=HBuf(J,I)/DSQRT(ATMASS(1+(J-1)/3))
633      CONTINUE
        ELSE
         CALL ZERO(HBUF(1,I),NRX)
        ENDIF
 623   CONTINUE
C
C Note that the REALATOM record is also written eleswhere. The normal
C modes in NORMCORD record is without mass weighing. The SYM_NMDS
C record essentially have the same info with out mass weights. May
C be we can use this record to store normal modes in internal coordiantes.
C Ajith Perera, 02/09.
       CALL IPUTREC(20,'JOBARC','REALATOM',1,INT(NRX/3))
       CALL DPUTREC(20,'JOBARC','NORMCORD',NRX*NRX,HBUF)
C
       IF(IFLAGS(68).NE.1) THEN
C
        CALL DGETREC(20,'JOBARC','BMATRIXC',NRX*NXM6,BMAT)
        CALL XGEMM('N','N',NXM6,NRX,NRX,1.0D0,BMAT,NXM6,HBUF,NRX,
     &             0.0D0,BUF1,NXM6)
        DO 624 I=1,NRX
         IF(TYPE(I).EQ.' VIBRATION ')THEN
          CALL NORMAL(BUF1(1,I),NXM6)
         ENDIF
624     CONTINUE
        IBOT=1
        WRITE(6,*)
        WRITE(6,3117)
3117    FORMAT(T20,'   Normal modes in internal coordinates ')
3114    ITOP=MIN(IBOT+4,NRX)
        WRITE(6,3116)
3116    FORMAT(75('-'))
        WRITE(6,'(5X,5(5X,F9.3))')(H(J,J),J=IBOT,ITOP)
        WRITE(6,3116)
        DO 3115 ICOORD=1,NXM6
#ifdef _DEBUG_LVLM1
         Print*, "In vib1: printing Normal Modes"
         Print*, VARNAM(ISQUASH(ICOORD))
         Print*, (ISQUASH(I), I=1,NXM6)
#endif
         WRITE(6,'(A,5(5X,F9.6))')VARNAM(ISQUASH(ICOORD)),
     &         (BUF1(ICOORD,J),J=IBOT,ITOP)
3115    CONTINUE
        IBOT=ITOP+1
        IF(IBOT.LT.NRX)GOTO 3114
        WRITE(6,3116)
        WRITE(6,*)
       ENDIF
C
C IF DIPDERS WERE EVALUATED, WRITE OUT DIPOLE MOMENT FUNCTION.
C
      IF(DIPOLE)THEN
       OPEN(UNIT=95,FILE='DIPSCR',FORM='UNFORMATTED',STATUS='OLD')
       READ(95)SCRATCH
       CLOSE(UNIT=95,STATUS='DELETE')
       CLOSE(95)
       WRITE(LUOUT,1111)
       WRITE(LUOUT,1109)
1109   FORMAT(T21,'Dipole Moment Function',/,T20,'(Normal ',
     & 'Coordinate Basis)')
       WRITE(LUOUT,1111)
       WRITE(LUOUT,1112)
1112   FORMAT(T5,'Mode',T14,'Symmetry',T26,'d(Mu(x))/dQ',T40,
     &'d(Mu(y))/dQ',T54,'d(Mu(z))/dQ')
       WRITE(LUOUT,1111)
       DO 1115 I=1,NRX
        IF(TYPE(I).EQ.' VIBRATION ')THEN
         IF(I/10.EQ.0)THEN
          WRITE(LUOUT,1121)I,IRREP(I),(SCRATCH(NRX*(J-1)+I),J=1,3)
         ELSEIF(I/100.EQ.0)THEN
          WRITE(LUOUT,1122)I,IRREP(I),(SCRATCH(NRX*(J-1)+I),J=1,3)
         ELSEIF(I/1000.EQ.0)THEN
          WRITE(LUOUT,1123)I,IRREP(I),(SCRATCH(NRX*(J-1)+I),J=1,3)
         ENDIF
1121     FORMAT(T6,'Q',I1,T16,A4,T27,F9.6,T41,F9.6,T55,F9.6)
1122     FORMAT(T6,'Q',I2,T16,A4,T27,F9.6,T41,F9.6,T55,F9.6)
1123     FORMAT(T5,'Q',I3,T16,A4,T27,F9.6,T41,F9.6,T55,F9.6)
        ENDIF
1115   CONTINUE
       WRITE(LUOUT,1111)
      ENDIF
C
C WRITE INFORMATION TO FREQUENCY ARCHIVE FILE FOR SUBSEQUENT
C  PICKUP BY STAT. THERMO. PROGRAM
C
C
C FIRST FILL SCRATCH VECTOR WITH "VIBRATIONAL" FREQUENCIES.
C
      ZMASS=0.D0
      DO 1039 I=1,NRX/3
      ZMASS=ZMASS+ATMASS(I)
 1039 CONTINUE
      CALL ZERO(SCRATCH,NRX)
      IX=0
      NIMAG=0
      DO 2231 I=1,NRX
      IF(TYPE(I).EQ.' VIBRATION ')THEN
       IX=IX+1
       SCRATCH(IX)=H(I,I)
       IF(FLAG(I).EQ.'i')NIMAG=NIMAG+1
      ENDIF
 2231 CONTINUE
      OPEN(UNIT=79,FILE='FRQARC',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(79)
      IREC=0
 2233 READ(79,END=2234,ERR=2234)LABEL
      IREC=IREC+1
      GOTO 2233
 2234 REWIND(79)
      DO 2251 I=1,IREC
       READ(79,ERR=2239,END=2239)LABEL
 2251 CONTINUE
      GOTO 6666
 2239 WRITE(LUOUT,781)
 781  FORMAT(T3,'@VIB1-W, Mysterious problem with FRQARC file.  ',
     &/,T5,' Existing version will be deleted.')
      CLOSE(UNIT=79,STATUS='DELETE')
      OPEN(UNIT=79,FILE='FRQARC',FORM='UNFORMATTED',STATUS='NEW')
 6666 LABEL='***FRQC***'
C
C The integer*8 NRXD and NIMAGD is introduced to maintain the compatibility
C 64-bit integer ACESII/a2proc, 12/10, Ajith Perera. 
C
      NRXD   = NRX
      NIMAGD = NIMAG
      WRITE(79)LABEL,FPGRP,NRXD,NIMAGD,ZMASS,ZIX,ZIY,ZIZ,
     &(SCRATCH(J),J=1,IX)
      CLOSE(UNIT=79,STATUS='KEEP')
C
C COMPUTE ROTATIONALLY PROJECTED FREQUENCIES ACCORDING TO THE MILLER-HAN
C METHOD
C
      ITIM=0
      CALL ZERO(SCRATCH,3*NX)
C
C CHECK TO SEE IF THIS IS A LINEAR MOLECULE.
C
      ILINER=0
      IF(FPGRP(2:2).EQ.'X')ILINER=1
c      DO 10 I=1,NRX
c      IPOS=1+ITIM*NRX
c      IF(TYPE(I).EQ.'TRANSLATION')THEN
c       CALL VADD(SCRATCH(IPOS),SCRATCH(IPOS),HBUF(1,I),NRX,1)
c       ITIM=ITIM+1
c      ENDIF
c10    CONTINUE
C
C COMPUTE ROTATIONAL PART OF PROJECTOR.
C
      DO 10 I=1,NRX-2,3
      NAT=1+(I-1)/3
      ZZ=DSQRT(ATMASS(NAT))
      SCRATCH(I)=0.D0
      SCRATCH(I+1)=-Q(I+2)*ZZ
      SCRATCH(I+2)=Q(I+1)*ZZ
      SCRATCH(NRX+I)=Q(I+2)*ZZ
      SCRATCH(NRX+I+1)=0.D0
      SCRATCH(NRX+I+2)=-Q(I)*ZZ
      SCRATCH(2*NRX+I)=-Q(I+1)*ZZ
      SCRATCH(2*NRX+I+1)=Q(I)*ZZ
      SCRATCH(2*NRX+I+2)=0.D0
 10   CONTINUE
      CALL NORMAL(SCRATCH,NRX)
      CALL NORMAL(SCRATCH(NRX+1),NRX)
      CALL NORMAL(SCRATCH(2*NRX+1),NRX)
      ITIM=3
      IF(ILINER.EQ.1)ITIM=2
      CALL ZERO(HBUF,NRX*NRX)
      CALL OUTERP(HBUF,SCRATCH,NRX,ITIM,0)
C
C NOW DO TRANSLATIONAL PART OF PROJECTOR.
C
      CALL ZERO(SCRATCH,3*NRX)
      DO 11 I=1,NRX-2,3
      ZZ=DSQRT(ATMASS(1+(I-1)/3))
      SCRATCH(I)=ZZ
      SCRATCH(NRX+I+1)=ZZ
      SCRATCH(2*NRX+I+2)=ZZ
 11   CONTINUE
      CALL NORMAL(SCRATCH,NRX)
      CALL NORMAL(SCRATCH(NRX+1),NRX)
      CALL NORMAL(SCRATCH(2*NRX+1),NRX)
      CALL ZERO(H,NRX*NRX)
      ITIM=3
      CALL OUTERP(H,SCRATCH,NRX,ITIM,1)
      CALL VADD (HBUF,H,HBUF,NRX*NRX,-1.D0)
      REWIND(78)
      READ(78)((H(I,J),J=1,NRX),I=1,NRX)
      CALL MODMATMUL(H,H,HBUF,NRX,NRX,NRX,NRX,NRX,NRX)
      WRITE(78)((H(I,J),J=1,NRX),I=1,NRX)
      CALL ZERO(H,NRX*NRX)
      CALL MTRANSP(HBUF,H,NRX,NRX,NRX,NRX)
      BACKSPACE(78)
      READ(78)((HBUF(I,J),J=1,NRX),I=1,NRX)
      CALL MODMATMUL(H,H,HBUF,NRX,NRX,NRX,NRX,NRX,NRX)
      CALL ZERO(HBUF,NRX*NRX)
      CALL EIG(H,HBUF,NRX,NRX,0)
c   o gather the diagonal eigenvalues back into column 1 for dumping to JOBARC
c     (Obviously we cannot use this matrix afterward.)
      do i = 2, nrx
         h(i,1) = h(i,i)
      end do
      call dputrec(20,'JOBARC','FORCECON',NRX,h)
      WRITE(LUOUT,9901)
9901  FORMAT(T3,' Vibrational frequencies after rotational projection ',
     &' of ',/,T3,' Cartesian force constants: ')
      FACT=5.14048D03
      ZPE=0.0D0
      DO 9902 I=1,NRX
      IF(H(I,I).LT.0.D0)THEN
       H(I,I)=FACT*DSQRT(-H(I,I))
       FLAG(I)='i'
      ELSE
       H(I,I)=FACT*DSQRT(H(I,I))
       FLAG(I)=' '
       ZPE=ZPE+H(I,I)
      ENDIF
      WRITE(LUOUT,9903)I,H(I,I),FLAG(I)
9903  FORMAT(T10,I4,T25,F10.4,A1)
9902  CONTINUE
      WRITE(LUOUT,9904)HALF*ZPE*FACT01
9904  FORMAT(T3,'Zero-point vibrational energy = ',F9.4,' kcal/mole.')
      ISTOP=100
      CLOSE(UNIT=78,STATUS='DELETE')
      CLOSE(UNIT=56,STATUS='KEEP')
      CALL CHKECK(NRX,HBUF,ECKART)
      RETURN
      END





