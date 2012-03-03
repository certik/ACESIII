
C THIS ROUTINE DRIVES THE FORMATION OF THE FORCE CONSTANTS.

c INPUT
c integer NATOM  : the number of atoms
c integer NIRREP : the number of irreps
c integer IORDGP :
c char*4  TYPE   : (FULL|COMP) point group
c integer NDSCR  : the amount of double scratch at DSCR

c OUTPUT
c char*8  LABEL(NIRREP)     : (scr) symmetry label of each irrep
c integer ISYMIRR(3*NATOM)  : (scr)
c integer INVOP(3*NATOM)    : (scr)
c double  SYOP(3*NATOM)     : (scr)
c double  SYMHESS(3*NATOM)  : Symmetry  coordinate Hessian
c double  CARTHESS(3*NATOM) : Cartesian coordinate Hessian
c double  DIPDER(9*NATOM)   : Symmetry  coordinate dipole derivative
c double  DIPDFUL(9*NATOM)  : Cartesian coordinate dipole derivative
c double  POLDER(27*NATOM)  : Symmetry  coordinate polarizability derivative
c double  POLDFUL(27*NATOM) : Cartesian coordinate polarizability derivative
c double  DSCR(NDSCR)       : (scr) double scratch

c RECORDS
c get 'NUMPOINT'
c get 'ENGPOINT'
c get 'REFENERG'
c get 'GRDPOINT'
c get 'DIPPOINT'
c get 'POLPOINT'
c get 'OPERSREF'
c get 'NUMVIBRT'
c get 'INVPSMAT'
c get 'NPTIRREP'
c get TYPE//'SYQT'
c get TYPE//'DEGN'
c get TYPE//'LABL'
c put 'HESSIANM'

      SUBROUTINE SETFCM(NATOM,NIRREP,IORDGP,TYPE,
     &                  LABEL,ISYMIRR,INVOP,SYOP,
     &                  SYMHESS,CARTHESS,
     &                  DIPDER,DIPDFUL,
     &                  POLDER,POLDFUL,
     &                  DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 TYPE
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ISYMIRR(3*NATOM),INVOP(3*NATOM),SYOP(9*IORDGP)
      DIMENSION SYMHESS(9*NATOM*NATOM),CARTHESS(9*NATOM*NATOM),
     &          DIPDER(3*3*NATOM),DIPDFUL(3*3*NATOM),
     &          POLDER(3*NATOM*9),POLDFUL(9*3*NATOM)
      double precision dscr(ndscr)

      DIMENSION idegen(100),nptsirr(100)
      LOGICAL PRINTQ, bTmp

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
#include "control.com"

      PRINTQ=(IFLAGS(1).GT.10)
      STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
      NSIZE=3*NATOM
      lFREE=1

      CALL GETREC(20,'JOBARC','NUMPOINT',1,NPOINT)
      IF (ENERONLY) THEN
         lENGPNT = lFREE
         lFREE   = lFREE + NPOINT
      ELSE
         lGRDPNT = lFREE
         lFREE   = lFREE + NSIZE*NPOINT
         lDIPPNT = lFREE
         lFREE   = lFREE + 3*NPOINT
         lPOLPNT = lFREE
         lFREE   = lFREE + 9*NPOINT
      END IF
      NDSCRLFT = NDSCR+1-lFREE
      IF (NDSCRLFT.LT.0) THEN
         print *, '@SETFCM: Insufficient memory.'
         print *, '         need ',-ndscrlft*ifltln,' more bytes'
         call aces_exit(1)
      END IF
      IF (ENERONLY) THEN
         CALL GETREC(20,'JOBARC','ENGPOINT',NPOINT*IINTFP,DSCR(lENGPNT))
         Print*, "Energy point"
         Print*, (DSCR(I), I=1, NPOINT)

#ifdef _ASSERT
         do i = 0, npoint-1
            if (dscr(lENGPNT+i).eq.0.d0) then
c            o if any energy is exactly 0, then ACES did not do all points
               print *, '@SETFCM: Assertion failed.'
               print *, '         Energy of point ',1+i,' is 0. a.u.'
               call aces_exit(1)
            end if
         end do
#endif
         CALL GETREC(20,'JOBARC','REFENERG',IINTFP,E0)
         Print*, E0
      ELSE
         CALL GETREC(20,'JOBARC','GRDPOINT',
     &               NSIZE*NPOINT*IINTFP,DSCR(lGRDPNT))
#ifdef _ASSERT
         do i = 0, nsize*npoint-1, nsize
            bTmp = .true.
            do j = 0, nsize-1
               if (dscr(lGRDPNT+i+j).ne.0.d0) bTmp = .false.
            end do
            if (bTmp) then
c            o if any gradient is exactly 0, then something is PROBABLY wrong
               print *, '@SETFCM: Assertion failed.'
               print *, '         Gradient of point ',1+(i/nsize),
     &                  ' is 0. a.u.'
               call aces_exit(1)
            end if
         end do
#endif
         CALL GETREC(20,'JOBARC','DIPPOINT',
     &               3*NPOINT*IINTFP,DSCR(lDIPPNT))
         CALL GETREC(20,'JOBARC','POLPOINT',
     &               9*NPOINT*IINTFP,DSCR(lPOLPNT))
         CALL GETREC(20,'JOBARC','OPERSREF',IORDGP*9*IINTFP,SYOP)
      END IF

      CALL GETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      CALL GETREC(20,'JOBARC','INVPSMAT',NMODE,INVOP)
      CALL GETREC(20,'JOBARC','NPTIRREP',NIRREP,NPTSIRR)
      CALL GETREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISYMIRR)
      CALL GETREC(20,'JOBARC',TYPE//'DEGN',NIRREP,IDEGEN)
      IF (PRINTQ) THEN
      CALL GETREC(20,'JOBARC',TYPE//'LABL',NIRREP*IINTFP,LABEL)
      END IF

C-----------------------------------------------------------------------
CJDW 6/21/96. Initialize other arguments. Try to make sure intensities
C             come out zero in ENERONLY jobs.
CJDW/AP Jul/98 extended to calculate Raman Intensities.
C-----------------------------------------------------------------------
      CALL ZERO(SYMHESS,9*NATOM*NATOM)
      CALL ZERO(DIPDER, 9*NATOM      )
      CALL ZERO(POLDER, 9*NATOM*3    )
C-----------------------------------------------------------------------

      IFIRST=0
      IFIRSTP=0
      IFIRSTD=0
      IFIRSTG=0
      IOFF=1
      IOFFD=1
      IOFFP=1
      IOFFINV=1
      NLEFT=NMODE

c   o loop over symmetry blocks
      DO IRREP=1,NIRREP

c   o find first occurance of this irrep
      ILOC=ISRCHEQ(NMODE,ISYMIRR,1,IRREP)
      IF (ILOC.NE.NMODE+1) THEN
         ILAST=ISRCHNE(NLEFT,ISYMIRR(ILOC),1,IRREP)
         NVIBSYM=ILAST-1
         NVIBUNQ=NVIBSYM/IDEGEN(IRREP)
         IF (PRINTQ) THEN
            WRITE(6,2000)LABEL(IRREP),IDEGEN(IRREP),NVIBSYM
2000        FORMAT(T3,' Symmetry : ',A,' Degeneracy : ',I1,
     &             ' Unique symmetry coordinates : ',I3)
         END IF

C USE FIRST NVIBUNQ VECTORS FOR DEGENERATE REPS SINCE THEY ARE SORTED
C ACCORDING TO SUBGROUP IRREPS.  OTHERS ARE REDUNDANT.

         IF (ENERONLY) THEN
#ifdef _DEBUG_LVL0
            Print*, "Entering ener2fcm"
#endif
            CALL ENER2FCM(NVIBUNQ,DSCR(lENGPNT+IFIRST),SYMHESS(IOFF),
     &                    INVOP(IOFFINV),STPSIZ,E0,DSCR(lFREE),NDSCRLFT)
            IFIRST=IFIRST+NPTSIRR(IRREP)
         ELSE
C 08/16, Extensions to do Raman Intensities, Ajith and John.
#ifdef _DEBUG_LVL0
            Print*, "Entering the grad2fcm"
#endif
            CALL GRAD2FCM(NATOM,NVIBUNQ,
     &                    DSCR(lGRDPNT+IFIRSTG),SYMHESS(IOFF),
     &                    DSCR(lDIPPNT+IFIRSTD),DIPDER(IOFFD),
     &                    DSCR(lPOLPNT+IFIRSTP),POLDER(IOFFP),
     &                    INVOP(IOFFINV),SYOP,STPSIZ,
     &                    DSCR(lFREE),NDSCRLFT)
            IOFFD=IOFFD+3*NVIBUNQ
            IOFFP=IOFFP+9*NVIBUNQ
            IFIRSTG=IFIRSTG+NPTSIRR(IRREP)*NSIZE+NVIBSYM
            IFIRSTD=IFIRSTD+NPTSIRR(IRREP)*3
            IFIRSTP=IFIRSTP+NPTSIRR(IRREP)*9
         END IF
         IOFF=IOFF+NVIBUNQ*NVIBUNQ
         IOFFINV=IOFFINV+NVIBSYM
         NLEFT=NLEFT-NVIBUNQ

      END IF
      END DO

c   o reset the scratch (l*PNT is not used anymore)
      lFREE = 1
      NDSCRLFT = NDSCR

c   o transform and write the Hessian to JOBARC
      CALL TRNFCM(NMODE,NATOM,NIRREP,
     &            SYMHESS,CARTHESS,
     &            ISYMIRR,IDEGEN,TYPE,DSCR(lFREE),NDSCRLFT)
      CALL PUTREC(20,'JOBARC','HESSIANM',NSIZE*NSIZE*IINTFP,CARTHESS)

C Modifications for Raman Intensities Ajith and John 08/98
      CALL TRNDIP(NMODE,NATOM,NIRREP,
     &            DIPDER,DIPDFUL,POLDER,POLDFUL,
     &            ISYMIRR,IDEGEN,TYPE,DSCR(lFREE),NDSCRLFT)

      RETURN
      END

