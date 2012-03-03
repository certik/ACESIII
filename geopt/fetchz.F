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

      SUBROUTINE FETCHZ(VERY_FIRST_TIME, Z, MAXMEM)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "linelen.par"
#include "coord.com"

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
      CHARACTER*(linelen) ZLINE, DUMSTR
      integer izline(80), izl(2,7)
      INTEGER ATOI, TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH 
      CHARACTER*10 VARTMP(3*MXATMS)
      CHARACTER*5 CHTEST
      CHARACTER*4 INTLABEL(MAXREDUNCO)
      CHARACTER*1 achar, czTmp, czPercent, czAsterisk
      Character*4 FPGrp, BPGrp, PGrp
      LOGICAL YESNO,BAD123,XYZIN,NWFINDIF,I_UNIQUE(MAXREDUNCO),
     &        VERY_FIRST_TIME, UNCONSTRAINED
      logical bStruct, Write_as_read
C
      double precision rtmp(3*mxatms)

C The arrays in FETCHZ that need to be managed dynamically are:
C
C REDUNCO: Store the RICs. It is of length 3*MXATMS.
C
C IBNDTO : Bonding connectivities for RICs. It is of length MXATMS*MXATMS.
C IREDUNCO : The connectivities for RICs. It is of length
C            4*MAXREDUNCO.
C EIGVECTORS: Keep the eigenvectors of the B-matrix. Needs to
C             to create the Gen. inverse of B-matrix. It is
C             of length MAXREDUNCO*MAXREDUNCO.
C
      DIMENSION JAX(3*MXATMS), JLIST(3*MXATMS), Z(MAXMEM),
     &          REDUNCO(3*MXATMS), IBNDTO(MXATMS*MXATMS), 
     &          IREDUNCO(4, MAXREDUNCO)
      DIMENSION EIGVECTORS(MAXREDUNCO, MAXREDUNCO)
C
#include "cbchar.com"
c
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     &                ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,
     &                IBASIS,XYZTol
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
#include "ric_heap.com"
C
      DATA BTOA /0.529177249D0/
      DATA ZILCH /0.0/
C
      ISQ(N) = MAX( N-6 , (N/3)+MOD(N,3)-1 )

      czPercent  = achar(_IACHAR_PERCENT)
      czAsterisk = achar(_IACHAR_ASTERISK)

c ----------------------------------------------------------------------

c   o fill COORD with zeroes
      if (VERY_FIRST_TIME) then
         call zero(q,3*MXATMS)
C         call zero(r, MAXREDUNCO)
         call zero(r, 3*MXATMS)
         call zero(atmass,MXATMS)
      end if

C     Symbol to indicate optimized coordinate
      CALL IGETREC(1,'JOBARC','FIRSTRUN',1,iTmp)
      PRINTI = (iTmp.NE.0)
      PRINTF = .TRUE.
      IONE = 1

#ifdef _DEBUG_LVL0
         Print*, "@-Entry to fetchz logic", 
     &           iFlags2(h_IFLAGS2_geom_opt)
#endif
C
      IF (IFLAGS(68).EQ.1) THEN
C
C The logic here pertains to the Cartesian Coordinate input. Until, 
C 2002, ACES was only capable of doing single point calculations
C with Cartesian coordinates. The internal redundent coordinate
C geometry optimizations with Cartesian coordinate input was
C long anticipated and much desired, but remanind only as such
C until 2000. During early 2000 Ajith Perera developed the routines
C that are needed for IRC, but they were remained not debugged or
C incorporated in to the general structure, since he had other 
C important commitments and had no help until summer 2002, when
C he finally got some help from Luis Galiano to complete what 
C started in 2000.
C  
         XYZIN = .TRUE.
         IF (iFlags2(h_IFLAGS2_geom_opt).eq.0) THEN
            Write(6,*) "Entering getxyz"
            CALL GETXYZ
            GO TO 3431
         END IF 
C          
         IF (iFlags2(h_IFLAGS2_geom_opt).eq.2) THEN
#ifdef _DEBUG_LVL0
         Print*, "Entering Cartesian Optimization block"
#endif
            CALL GETXYZ                     
C
#ifdef _DEBUG_LVL0
            Write(6, "(a,I5)") "The # of atoms", Natoms
            ioff = 0
            Do I = 1, Natoms
               Write(6, "(3F10.5)") (Q(ioff+j), J=1,3)
               ioff = ioff + 3
            Enddo
#endif
 
            CALL DCOPY(3*NATOMS, Q, 1, R, 1)
            ISCRATCH = 1
            ITMPQ    = ISCRATCH + 9*NATOMS
            INEWQ    = ITMPQ    + 3*NATOMS
            ILAST    = INEWQ    + 3*NATOMS
            CALL SYMMETRY(Z(ISCRATCH), Z(ITMPQ), Z(INEWQ), .TRUE.)
C
#ifdef _DEBUG_LVL0
               Write(6, "(3F10.5)") (R(j), J=1,3*NATOMS)
#endif
            NXM6 = 3*NATOMS 
            NX   = 3*NATOMS
C
c            CALL IZERO(NCON, 3*MXATMS)
            do i = 1, 3*MXATMS
               NCON(i) = 0
            enddo

            CALL ASSIGN_XYZ_LABELS(3*NATOMS, INTLABEL)
C
            NUNIQUE = 0
            NOPT    = 0
C
            DO IVAR = 1, 3*NATOMS
               I_UNIQUE(IVAR) = .TRUE.
               IF (I_UNIQUE(IVAR)) THEN
                   NUNIQUE = NUNIQUE + 1
                   IUNIQUE(NUNIQUE) = IVAR
               END IF
C
C This is the place where info about constrained optimizations
C are added. At the moment all the unique coordinates are
C optimized.   

               UNCONSTRAINED  = .TRUE.
               IF (UNCONSTRAINED) THEN
                  IF (I_UNIQUE(IVAR)) THEN
                      NOPT         = NOPT  + 1
                      VARNAM(IVAR) = INTLABEL(IVAR)
                      PARNAM(NOPT) = INTLABEL(IVAR)
                      NOPTI(NOPT)  = IVAR
                  ELSE
                      VARNAM(IVAR) = INTLABEL(IVAR)
                  ENDIF
               ENDIF
            END DO
C
#ifdef _DEBUG_LVL0
            WRITE(6,"(a,3(1x,I3))") "The number of Degrees of Fre =", 
     &                  NOPT, NUNIQUE, 3*NATOMS
            Write(6,*)
            Write(6, "(a4)") (Varnam(I), I=1,NOPT)
#endif
            DO I = 1, 3*NATOMS
               ISQUASH(I) = I
            END DO
C
C We have to rebuilt the NEQ array so that it represents the number
C of symmetry equivalent atoms for given symmetry unique atoms. The
C NEQ that comes out from PROCESS_REDINTRNLS corresponds to the
C natural order (ordered according to the unique degree of freedoms
C regardless of the atoms involved)
C
            DO ICOORD = 1, 3*NATOMS
               NCON(IUNIQUE(ICOORD)) = NEQ(ICOORD)
            ENDDO
C
c            CALL IZERO(NEQ, 3*NATOMS)
            do i = 1, 3*natoms
               NEQ(i) = 0
            enddo 
            CALL ICOPY(3*NATOMS, NCON, 1, NEQ, 1)
C
#ifdef _DEBUG_LVL0
         WRITE(6,*) "@FETCHZ the NEQ"
         WRITE(6,"(5I2)") (NCON(IIIK), IIIK = 1, 3*NATOMS)
     &        
#endif
            GO TO 3431
C
         END IF      
C
         IF (iFlags2(h_IFLAGS2_geom_opt).ge.3) THEN
C
C Redundent Internal Optimizations Starts Here
C
c      o allocate the RIC heap and create the pointers

c----------------------------------------------------------------------------
c   INIT_RIC_HEAP has been modified to read the JOBARC for REDNCORD to 
c   determine the actual memory allocation.  If VERY_FIRST_TIME is .TRUE.,
c   this info is not available.  So we must initialize it with the maxredunco
c   value in this case.
c
c   Successive calls to INIT_RIC_HEAP will read the actual number of 
c   redundant coordinates from JOBARC.
c----------------------------------------------------------------------------

         IF (VERY_FIRST_TIME) THEN
            TOTREDNCO = maxredunco
            CALL IPUTREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
         ENDIF

         IF (VERY_FIRST_TIME) THEN
            CALL GETXYZ
            call init_ric_heap
            ISCRATCH = 1
            ITMPQ    = ISCRATCH + 9*NATOMS
            INEWQ    = ITMPQ    + 3*NATOMS
            ILAST    = INEWQ    + 3*NATOMS
            CALL SYMMETRY(Z(ISCRATCH), Z(ITMPQ), Z(INEWQ), .TRUE.)
#ifdef _DEBUG_LVL0
            Print*, "Entering the Redundent Internal block"
#endif
            CALL GEN_REDNINTRNLS(Q, REDUNCO, INTLABEL, IBNDTO,
     &                           dRICHeap(z_BMat), dRICHeap(z_GMat),
     &                           TOTREDNCO, IATNUM, NATOMS, IREDUNCO,
     &                           TOTNOFBND, TOTNOFANG, TOTNOFDIH,
     &                           dRICHeap(z_BTGInv), EIGVECTORS)
C
            CALL IPUTREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
            CALL IPUTREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO, 
     &                  IREDUNCO)
            CALL IPUTREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
            CALL IPUTREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
            CALL IPUTREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
            CALL IPUTREC(20, 'JOBARC', 'IBONDTO ', NATOMS*NATOMS, 
     &                  IBNDTO)
C
            IF (IFLAGS(h_IFLAGS_symmetry).EQ. 2 .AND. 
     &          FPGRP .NE. 'C1') THEN
               CALL PROCESS_REDINTRNLS(Q, REDUNCO, IREDUNCO, IBNDTO,
     &                                 I_UNIQUE, FPGRP, TOTNOFBND,
     &                                 TOTNOFANG, TOTNOFDIH, TOTREDNCO,
     &                                 NEQ, IEQUIV, NATOMS, NCON, IPRNT,
     &                                 EPSILON)
               CALL IPUTREC(20,'JOBARC','PLSMINSP',TOTREDNCO,NCON)
            ELSE
c               CALL IZERO(NCON,TOTREDNCO)
c               CALL IZERO(NEQ, TOTREDNCO)
              
               DO IREDNCO = 1, TOTREDNCO
                  I_UNIQUE(IREDNCO) = .TRUE.
                  NCON(IREDNCO)     = 0
                  NEQ(IREDNCO)      = 0 
               ENDDO
            ENDIF
            CALL IPUTREC(20,'JOBARC','PLSMINSP',TOTREDNCO,NCON) 
            CALL IPUTREC(20, 'JOBARC', 'UNIQUEDF', TOTREDNCO, 
     &                  I_UNIQUE)
         ELSE
C
            CALL IGETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
            CALL IGETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO, 
     &                  IREDUNCO)
            CALL IGETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
            CALL IGETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
            CALL IGETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
C
            CALL BUILT_BGMTRX(Q, REDUNCO, IREDUNCO, TOTREDNCO,
     &                        TOTNOFBND, TOTNOFANG, TOTNOFDIH,
     &                        NATOMS,dRICHeap(z_BMat),dRICHeap(z_GMat),
     &                        EIGVECTORS,EPSILON,dRICHeap(z_BTGInv),
     &                        dRICHeap(z_DerBMat))
         ENDIF
C
C At this point, we need to update the /COORD/, /CBCHR/ common blocks  
C so that a geometry optimization can proceed with using the standard
C code. 
C  
         LENGTH_BGMAT = TOTREDNCO*NATOMS*3
         CALL DPUTREC(20,'JOBARC', 'BMATRIXT', LENGTH_BGMAT,
     &               dRICHeap(z_BMat))
         CALL DPUTREC(20,'JOBARC', 'GMATRIX ', LENGTH_BGMAT,
     &               dRICHeap(z_GMat))
C
C Copy the intial values of the RIC to R array. Before that we need 
C mark the plus-minus pairs.
C        
         IF (IFLAGS(h_IFLAGS_symmetry).EQ. 2) THEN
            CALL IGETREC(20,'JOBARC','PLSMINSP',TOTREDNCO,NCON)
            DO ICOORD = 1, TOTREDNCO
               IF (NCON(ICOORD).EQ.1) REDUNCO(ICOORD) = -REDUNCO(ICOORD)
            END DO
         END IF
C
         CALL DCOPY(TOTREDNCO, REDUNCO, 1, R, 1)
         NXM6 = TOTREDNCO
         NX   = MAX(TOTREDNCO, 3*NATOMS)
C
         IF (.NOT. VERY_FIRST_TIME) RETURN
C
c         CALL IZERO(NCON, 3*MXATMS)
         do i = 1, 3*mxatms
            ncon(i) = 0
         enddo
C
         NUNIQUE = 0
         NOPT    = 0
C
         DO IVAR = 1, TOTREDNCO
C
            IF (I_UNIQUE(IVAR)) THEN
                NUNIQUE = NUNIQUE + 1
                IUNIQUE(NUNIQUE) = IVAR
            END IF
C
C This is the place where info about constrained optimizations
C are added. At the moment all the unique coordinates are 
C optimized.
C 
            UNCONSTRAINED  = .TRUE.  
            IF (UNCONSTRAINED) THEN
               IF (I_UNIQUE(IVAR)) THEN
                  NOPT         = NOPT  + 1
                  VARNAM(IVAR) = INTLABEL(IVAR)
                  PARNAM(NOPT) = INTLABEL(IVAR)
                  NOPTI(NOPT)  = IVAR 
               ELSE
                  VARNAM(IVAR) = INTLABEL(IVAR) 
               ENDIF
            ENDIF
C
         END DO 
C
#ifdef _DEBUG_LVL0
         WRITE(6,"(a,3(1x,I3))") "The number of Degrees of Fre =", 
     &         NOPT, NUNIQUE, TOTREDNCO
#endif
C
         DO I = 1, TOTREDNCO
            ISQUASH(I) = I
         END DO
C 
C We have to rebuilt the NEQ array so that it represents the number
C of symmetry equivalent atoms for given symmetry unique atoms. The
C NEQ that comes out from PROCESS_REDINTRNLS corresponds to the 
C natural order (ordered according to the unique degree of freedoms
C regardless of the atoms involved)
C
         DO ICOORD = 1, TOTREDNCO
            NCON(IUNIQUE(ICOORD)) = NEQ(ICOORD)
         ENDDO
C
c         CALL IZERO(NEQ, TOTREDNCO)
         do i = 1, totrednco
            neq(i) = 0
         enddo  
         CALL ICOPY(TOTREDNCO, NCON, 1, NEQ, 1)
C
#ifdef _DEBUG_LVL0
         WRITE(6,*) "@FETCHZ the NEQ"
         WRITE(6,"(5I3,I)") (NCON(IIIK), IIIK = 1, TOTREDNCO),
     &         TOTREDNCO
#endif
         NXM6 = TOTREDNCO
         NX   = MAX(TOTREDNCO, 3*NATOMS)
C 
         GO TO 3431
C
C Endif for GEOM_OPT=RIC
       ENDIF
C
C Endif for COORD=CARTESIAN
      END IF
C
c      CALL IZERO(JLIST,3*MXATMS)
      do i = 1, 3*mxatms
         jlist(i) = 0
         jax(i)   = 0 
      enddo 
c      CALL IZERO(JAX,3*MXATMS)
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
            write(*,*) '@FETCHZ: The number of atoms has exceeded the ',
     &                 'maximum allowed of ',mxatms
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
     &         '@FETCHZ: The Z-matrix was not terminated properly.'
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
         IF (PRINTI) WRITE(LUOUT,*) '@FETCHZ: The ordering of the ',
     &      'first two atoms in the Z-matrix has been changed.'
      END IF

      IF (PRINTI .AND. IPRNT.GE.2) THEN
         WRITE (LUOUT,*) '@FETCHZ: plain Z-matrix'
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
c      CALL IZERO(NEQ,3*MXATMS)
      do i = 1, 3*mxatms
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
      call rd_intcoord(R)

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
            call dgetrec(0,'JOBARC','RXSTRUCT',J,rtmp)
            if (j.gt.0) then
               call dgetrec(1,'JOBARC','RXSTRUCT',J,rtmp)
               write(luout,1002)
               DO I = 1, NUNIQUE
                  WRITE(LUOUT,'(T25,A,T38,F11.6)')
     &               VARNAM(ISQUASH(IUNIQUE(I))),
     &               RTMP(ISQUASH(IUNIQUE(I)))
               END DO
            end if
c         o product
            call dgetrec(0,'JOBARC','PRSTRUCT',J,rtmp)
            if (j.gt.0) then
               call dgetrec(1,'JOBARC','PRSTRUCT',J,rtmp)
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
C
C GET ATOMIC MASSES, ATOMIC NUMBERS FROM PERIODIC TABLE SUBROUTINE
C
      CALL IPUTREC(20, 'JOBARC', 'CONCTVTY', 3*NATOMS, NCON)
      CALL IPUTREC(20, 'JOBARC', 'ICSQUASH', 3*NATOMS, ISQUASH)
      CALL PUTCREC(20, 'JOBARC', 'INTCNAM', 15*NATOMS, VARNAM)

      Call IGETREC(20,  'JOBARC', 'HAVEGEOM', 1, I_havegeom)
      Call IGETREC(0, 'JOBARC', 'PASS1   ', Length, Itmp)
      Write_as_read = (Length .LT. 0)
#ifdef _DEBUG_LVL0
      Print*, "Write_as_read from FETCHZ:", Write_as_read 
      Print*, "Internal coordinates in FETCHZ"
      Print*,  (R(I), I=1, 3*NATOMS)
#endif
      Call Angs2bohr(R, 3*NATOMS)
      If (Write_as_read) CALL DPUTREC(20, 'JOBARC', 'CORD_INT', 
     &                               3*NATOMS, R)
      Call Bohr2Angs(R, 3*NATOMS)  
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
         CALL DGETREC(20,'JOBARC','NEXTGEOM',3*NATOMS,Q)
      END IF
      Write(6, "(3F10.5)")(Q(I), I=1, 3*NATOMS)

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

      RETURN
      END

