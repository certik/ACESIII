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

C Read the Z-matrix and associated variables

#include "iachar.h"

      SUBROUTINE read_zmat_internal_coords()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "linelen.par"

      LOGICAL PRINTI,PRINTF,BTMP

c     This number is supposed to be a unique tag to see whether some
c     internal coordinates were not defined.
      DOUBLE PRECISION dSmallNumber
      PARAMETER (dSmallNumber=5.361d-21, EPSILON= 1.0D-07)

#include "io_units.par"
#include "flags.h"

C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      CHARACTER*(linelen) ZLINE, DUMSTR
      integer izline(80), izl(2,7)
      INTEGER ATOI, TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH 
      CHARACTER*10 VARTMP(3*MXATMS)
      CHARACTER*5 CHTEST
      CHARACTER*4 INTLABEL(3*MXATMS)
      CHARACTER*1 achar, czTmp, czPercent, czAsterisk
      Character*4 FPGrp, BPGrp, PGrp
      LOGICAL YESNO,BAD123,XYZIN,NWFINDIF,I_UNIQUE(MAXREDUNCO),
     &        VERY_FIRST_TIME, UNCONSTRAINED
      logical bStruct
C
      double precision rtmp(3*MXATMS)

C
      DIMENSION JAX(3*MXATMS), JLIST(3*MXATMS)
C
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(3*MXATMS),
     &     PARNAM(3*MXATMS)
#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     &                ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,
     &                IBASIS,XYZTol
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
C
      DATA BTOA /0.529177249D0/
      DATA ZILCH /0.0/
C
      ISQ(N) = MAX( N-6 , (N/3)+MOD(N,3)-1 )

      czPercent  = achar(_IACHAR_PERCENT)
      czAsterisk = achar(_IACHAR_ASTERISK)

c ----------------------------------------------------------------------

C     Symbol to indicate optimized coordinate
      CALL IGETREC(1,'JOBARC','FIRSTRUN',1,iTmp)
      PRINTI = (iTmp.NE.0)
      PRINTF = .TRUE.
      IONE = 1
C
      do i = 1, 3*MXATMS
         jlist(i) = 0
         jax(i)   = 0
         vartmp(i) = ' '
      enddo
C
      XYZIN   = .FALSE.
      BAD123  = .FALSE.
      IERR    = 0
      RTOA    = 180.0D0/DACOS(-1.0D0)
      NUNIQUE = 0
      NOPT    = 0

c ----------------------------------------------------------------------

c   o open up the main input file
      open(LuZ,FILE=ZFil,STATUS='OLD')
      rewind LuZ

c   o skip the header (zline=TITLE on end do)
      btmp = .true.
      do while (btmp)
         read(luz,'(a)') zline
         call parsez(zline,izl)
         i = izl(1,1)
         btmp = (i.eq.0).or.(zline(i:i).eq.czPercent)
      end do

c   o count the number of elements
      btmp = .true.
      natoms = 0
      do while (btmp)
         read(luz,'(a)') zline
         if (linblnk(zline).eq.0) then
            btmp = .false.
         else
            natoms = natoms + 1
         end if
         if (natoms.gt.mxatms) then
            write(*,*)
            write(*,*) '@READ_ZMAT: The number of atoms has ',
     &                 'exceeded the maximum allowed of ',mxatms
            write(*,*)
            call errex
         end if
      end do
      if (printi)
     &   write(*,'(a,i3,a)') ' ',natoms,' entries found in Z-matrix '

c     set number of coordinates based on number of elements found in Z-matrix
      NX = 3 * NATOMS

c Create an array that relates TR invariant (squashed) coords to full
c space for use as lookup during compression/decompression.
      ISQUASH(1) = 4
      ISQUASH(2) = 7
      ISQUASH(3) = 8

c Set NXM6 to 3*NATOMS-6 or just 1 if it is a diatomic.
      NXM6 = MAX(NX-6,1)
      DO I = 4,NXM6
         ISQUASH(I) = I+6
      END DO

c ----------------------------------------------------------------------

c Retrieve the Z-matrix:
c     atomic symbols into ZSYM
c     variable names into VARTMP
c     connectivity into NCON
c
c Do this by reading in each line as text and picking it apart using
c the positions of white space determined in PARSEZ.

c   o skip the header (zline=TITLE on end do)
      rewind(luz)
      btmp = .true.
      do while (btmp)
         read(luz,'(a)') zline
         call parsez(zline,izl)
         i = izl(1,1)
         btmp = (i.eq.0).or.(zline(i:i).eq.czPercent)
      end do

c   o save the title
      if (printi) then
         i = linblnk(zline)
         if (i.lt.68) then
            write(LuOut,*)
            write(LuOut,'(2a)') ' Job Title : ', zline(1:i)
            write(LuOut,*)
         else
            write(LuOut,*)
            write(LuOut,'(1a)') ' Job Title : '
            write(LuOut,'(a)') zline(1:i)
            write(LuOut,*)
         end if
      end if
c YAU : old
c      read(zline,'(80a1)') izline
c      call iputrec(1,' ','TITLE',80,izline)
c YAU : new
      call putcrec(1,' ','TITLE',80,zline(1:80))
c YAU : end

c   o first element
      read(luz,'(a)') zline
      call parsez(zline,izl)
      if (izl(1,1).eq.0.or.izl(1,2).ne.0) call debug_zmat
      ZSYM(1) = ZLINE(izl(1,1):izl(2,1))
      VARTMP(1) = ' '
      VARTMP(2) = ' '
      VARTMP(3) = ' '
      NCON(1)   = 0
      NCON(2)   = 0
      NCON(3)   = 0
      IF (NATOMS.EQ.1) THEN
         NXM6 = 0
c SB: Fixed output for atomic calculations
         nunique = 0
         nopt = 0
         goto 3441
      END IF

c   o second element
      read(luz,'(a)') zline
      call parsez(zline,izl)
      if (izl(1,3).eq.0.or.izl(1,4).ne.0) call debug_zmat
      ZSYM(2)   = ZLINE(izl(1,1):izl(2,1))
      VARTMP(4) = ZLINE(izl(1,3):izl(2,3))
      VARTMP(5) = ' '
      VARTMP(6) = ' '
      NCON(4)   = 1
      NCON(5)   = 0
      NCON(6)   = 0
      IF (NATOMS.EQ.2) GOTO 34

c   o third element
      read(luz,'(a)') zline
      call parsez(zline,izl)
      if (izl(1,5).eq.0.or.izl(1,6).ne.0) call debug_zmat
      ZSYM(3)   = ZLINE(izl(1,1):izl(2,1))
      VARTMP(7) = ZLINE(izl(1,3):izl(2,3))
      VARTMP(8) = ZLINE(izl(1,5):izl(2,5))
      VARTMP(9) = ' '
      NCON(7)   = 2
      NCON(8)   = 1
      NCON(9)   = 0

c SET A FLAG IF THIS IS NOT A 2--1--3 Z-MATRIX.
      JUNK   = ATOI( ZLINE(izl(1,2):izl(2,2)) )
      BAD123 = (JUNK.EQ.1)

c Ajith Perera 07/2000
c The extended (to handle MOLCAS integrals) vmol2ja needs to know whether
c the first two atoms have been switched.
      IF (BAD123) CALL IPUTREC(20,'JOBARC','12SWITCH',1,JUNK)

c   o remaining elements
      do i = 4, natoms
         read(luz,'(a)') zline
         if (zline(1:1).eq.czAsterisk) then
            write(*,*)
     &         '@READ_ZMAT: The Z-matrix was not terminated properly.'
            call errex
         end if
         call parsez(zline,izl)
         if (izl(1,7).eq.0) call debug_zmat
         ZSYM(I)       =       ZLINE(izl(1,1):izl(2,1))
         VARTMP(3*I-2) =       ZLINE(izl(1,3):izl(2,3))
         VARTMP(3*I-1) =       ZLINE(izl(1,5):izl(2,5))
         VARTMP(3*I)   =       ZLINE(izl(1,7):izl(2,7))
         NCON(3*I-2)   = ATOI( ZLINE(izl(1,2):izl(2,2)) )
         NCON(3*I-1)   = ATOI( ZLINE(izl(1,4):izl(2,4)) )
         NCON(3*I)     = ATOI( ZLINE(izl(1,6):izl(2,6)) )
c      o catch basic errors and immediately dump to debug_zmat()
         J = 3*I-2
         IF ( (NCON(J  ).EQ.NCON(J+1)) .OR.
     &        (NCON(J  ).EQ.NCON(J+2)) .OR.
     &        (NCON(J+1).EQ.NCON(J+2))      ) call debug_zmat
         IF ( (NCON(J  ).EQ.0) .OR.
     &        (NCON(J+1).EQ.0) .OR.
     &        (NCON(J+2).EQ.0)      ) call debug_zmat
         IF ( (NCON(J  ).GE.I) .OR.
     &        (NCON(J+1).GE.I) .OR.
     &        (NCON(J+2).GE.I)      ) call debug_zmat
c     end do i = 4, natoms
      end do

c IF CONNECTIVITY OF THE FIRST THREE ATOMS IS 3--1--2. SWITCH
c ATOMS 1 AND 2 AND THEN ALL REFERENCES TO THEM.
   34 CONTINUE
      IF (BAD123) THEN
         CHTEST  = ZSYM(1)
         ZSYM(1) = ZSYM(2)
         ZSYM(2) = CHTEST
         DO J = 9, NX
            ITMP = 0
            IF (NCON(J).EQ.2) ITMP = 1
            IF (NCON(J).EQ.1) ITMP = 2
            IF (ITMP.NE.0)    NCON(J) = ITMP
         END DO
         IF (PRINTI) WRITE(LUOUT,*) '@READ_ZMAT: The ordering of the ',
     &      'first two atoms in the Z-matrix has been changed.'
      END IF

      IF (PRINTI .AND. IPRNT.GE.2) THEN
         WRITE (LUOUT,*) '@READ_ZMAT: plain Z-matrix'
         J = 1
         DO I = 1, NAtoms
            WRITE (LUOUT,'(I2,A,1X,A,3(1X,I3,1X,A))')
     &            I,':',ZSYM(I),
     &            NCON(J+0),VARTMP(J+0),
     &            NCON(J+1),VARTMP(J+1),
     &            NCON(J+2),VARTMP(J+2)
            J = J + 3
         END DO
         WRITE (LUOUT,*)
      END IF

c REDUNDANCY CHECKING AND GENERATE VECTORS CONTAINING POSITIONS
c OF COORDINATES TO BE OPTIMIZED AND CORRESPONDING REDUNDANT
c COORDINATES. LISTS USE INTERNAL COORDINATE NUMBERING WHICH
c CORRESPONDS TO THE "SQUASHED" SCHEME (1,2,3,4,... vs 4,7,8 etc)

      do i = 1, 3*NATOMS
         neq(i) = 0
      enddo

      DO 210 I = 1, NX
         ISQ_I = ISQ(I)
         IF (.NOT.(I.LT.4.OR.I.EQ.5.OR.I.EQ.6.OR.I.EQ.9)) THEN
            IRED = 0
            DO J = 1, I-1
               IF (VARTMP(J).EQ.VARTMP(I)) THEN
                  IRED = 1
                  ISQ_J = ISQ(J)
                  NEQ(ISQ_J) = NEQ(ISQ_J) + 1
                  IEQUIV(JLIST(ISQ_J),NEQ(ISQ_J)) = ISQ_I
               END IF
            END DO
            IF (IRED.LE.0) THEN
               NUNIQUE = NUNIQUE + 1
               IUNIQUE(NUNIQUE) = ISQ_I
               JLIST(ISQ_I) = NUNIQUE
            END IF
C
C CHECKING FOR PARAMETERS TO BE OPTIMIZED
C
            JAX(I) = INDEX(VARTMP(I),czAsterisk)
            VARNAM(I) = VARTMP(I)
            IF (JAX(I).NE.0.OR.iFlags2(h_IFLAGS2_geom_opt).EQ.2) THEN
               IF (JAX(I).NE.0) VARNAM(I) = VARTMP(I)(:JAX(I)-1)
               IF (IRED.EQ.1) GOTO 210
               NOPT = NOPT + 1
               PARNAM(NOPT) = VARNAM(I)
               NOPTI(NOPT) = ISQ_I
            END IF
         END IF
  210 CONTINUE

c   o read the blank line after the Z-matrix
      READ(LUZ,'(A)') ZLINE

c   o load the first structure into R
      call rd_intcoord(r)
      print *,'R after RD_INTCOORD: '
      do i = 1, 3*natoms
         print *,'i, r(i) ',i,r(i)
      enddo

c   o read additional structures until we hit a namelist
c     (iStruct corresponds to the structure we expect to read)
      iStruct = 1
      bStruct = .true.
      do while (bStruct.and.iStruct.lt.3)
         iStruct = iStruct + 1
         read(luz,'(a)') zline
         backspace(luz)
         ISTP = INDEX(zline(1:),'=')
         I    = INDEX(zline(1:),'*')
         bStruct = ((i.eq.0).and.(istp.ne.0))
         if (bStruct) then
            call rd_intcoord(rtmp)
            if (iStruct.eq.2) then
c            o read the next line and see if this structure was the TS or PR
               read(luz,'(a)') zline
               backspace(luz)
               ISTP = INDEX(zline(1:),'=')
               I    = INDEX(zline(1:),'*')
               btmp = ((i.ne.0).or.(istp.eq.0))
            else
               btmp = .true.
            end if
            if (btmp) then
               call dputrec(1,'JOBARC','PRSTRUCT',nx,rtmp)
            else
c            o move R to RXSTRUCT since we are optimizing the TS structure
               call dputrec(1,'JOBARC','RXSTRUCT',nx,R)
               call dcopy(nx,rtmp,1,R,1)
            end if
c        end if (bStruct)
         end if
c     end do while (bStruct)
      end do

 3441 continue

C WRITE OUT RELEVANT INFORMATION
      IF (PRINTI.AND. ncycle.EQ.0) THEN
         WRITE(LUOUT,*)
     &      ' There are ',NUNIQUE,' unique internal coordinates.'
         WRITE(LUOUT,*)
     &      ' Of these, ',NOPT ,' will be optimized.'
      END IF
      IF (IPRNT.GE.2.AND.PRINTI) THEN
         DO I = 1, NOPT
            WRITE(LUOUT,*) ' Parameter #',I,' to be optimized is ',
     &         PARNAM(I)
         END DO
      END IF

C PRINT OUT Z-MATRIX SO THAT USER CAN RECONSTRUCT LATER IF NEED BE.
      IF (IPRNT.GE.0.AND.PRINTI) THEN
         WRITE(LUOUT,1001)
 1001    FORMAT(T3,' User supplied Z-matrix: ',/,80('-'),/
     &          T8,'SYMBOL',T18,'BOND',
     &          T28,'LENGTH',T38,'ANGLE',T48,'ANGLE',T58,'DIHED',T68,
     &              'ANGLE',/,
     &          T18,' TO ',T27,'(ANGST)',T37,' WRT ',T47,'(DEG)',T57,
     &              ' WRT ',
     &          T67,'(DEG)',/,80('-'))
         WRITE(LUOUT,734) ZSYM(1)
c SB fixed output for atomic and diatomic calculations
         if (natoms.gt.1) WRITE(LUOUT,734) ZSYM(2),NCON(4),VARNAM(4)
         if (natoms.gt.2) WRITE(LUOUT,734) ZSYM(3),NCON(7),VARNAM(7),
     &                                     NCON(8),VARNAM(8)
         DO I = 4, NATOMS
            IT3 = 3 * I
                          WRITE(LUOUT,734) ZSYM(I),
     &                       (NCON(IT3-J),VARNAM(IT3-J),J=2,0,-1)
         END DO
  734    FORMAT(T9,A,T19,I2,T30,A,T40,I2,T50,A,T60,I2,T70,A)
         if (nunique.gt.0) then
            WRITE(LUOUT,1000)
 1000       FORMAT(T18,' *Initial values for internal coordinates* ',/,
     &             T23,'Name',t40,'Value')
 1002       FORMAT(T18,' *Initial values for reactant coordinates* ',/,
     &             T23,'Name',t40,'Value')
 1003       FORMAT(T18,' *Initial values for product  coordinates* ',/,
     &             T23,'Name',t40,'Value')
            DO I = 1, NUNIQUE
               WRITE(LUOUT,'(T25,A,T38,F11.6)')
     &            VARNAM(ISQUASH(IUNIQUE(I))),R(ISQUASH(IUNIQUE(I)))
            END DO
c         o reactant
            call getrec(0,'JOBARC','RXSTRUCT',J,rtmp)
            if (j.gt.0) then
               call getrec(1,'JOBARC','RXSTRUCT',J,rtmp)
               write(luout,1002)
               DO I = 1, NUNIQUE
                  WRITE(LUOUT,'(T25,A,T38,F11.6)')
     &               VARNAM(ISQUASH(IUNIQUE(I))),
     &               RTMP(ISQUASH(IUNIQUE(I)))
               END DO
            end if
c         o product
            call getrec(0,'JOBARC','PRSTRUCT',J,rtmp)
            if (j.gt.0) then
               call getrec(1,'JOBARC','PRSTRUCT',J,rtmp)
               write(luout,1003)
               DO I = 1, NUNIQUE
                  WRITE(LUOUT,'(T25,A,T38,F11.6)')
     &               VARNAM(ISQUASH(IUNIQUE(I))),
     &               RTMP(ISQUASH(IUNIQUE(I)))
               END DO
            end if
         END IF
         WRITE(LUOUT,735)
 735  FORMAT(80('-'))
      END IF

c ----------------------------------------------------------------------

cSSS 3431 continue

C GET ATOMIC MASSES, ATOMIC NUMBERS FROM PERIODIC TABLE SUBROUTINE
C
      CALL PERTABLE
C
 3431 continue
C
C FOR FINDIF CALCULATIONS, READ IN NEW CARTESIAN COORDINATES AND
C TRANSFORM THESE BACK INTO MODIFIED INTERNALS
C 
      IONE = 1
      CALL IGETREC(-1,'JOBARC','PASS1   ',IONE,INWFND)
      IF (INWFND.EQ.1) THEN
         XYZIN = .TRUE.
         IFLAGS(68) = 1
         CALL GETREC(20,'JOBARC','NEXTGEOM',3*NATOMS*IINTFP,Q)
      END IF

C COUNT NUMBER OF REAL ATOMS
C Also add a record that contains the sum of both real and dummy atoms.
C This is needed for external programs that use ACES II JOBARC records
C most of which are based on the total number of atoms, not just the
C number of real atoms. 01/2006, Ajith Perera.
C
      CALL IPUTREC(20,'JOBARC','ZMATATMS',IONE,NATOMS)
      NREAL = 0
      DO IATOM = 1, NATOMS
         IF (ATMASS(IATOM).NE.ZILCH) NREAL = NREAL + 1
      END DO
      CALL IPUTREC(20,'JOBARC','NREALATM',IONE,NREAL)

c   
c   Convert to Cartesian coords.
c
      call gmetry()
      print *,'Cartesian coordinates after gen_cart_coord:'
      do i = 1, 3*natoms
         print *,'i, q(i) ',i,q(i)
      enddo

      RETURN
      END

