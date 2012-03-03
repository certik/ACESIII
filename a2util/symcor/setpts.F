
C THIS ROUTINE DETERMINES WHICH POINTS WILL BE RUN AND WHAT TYPE
C OF CALCULATION (GRADIENT OR ENERGY ONLY) IS MOST EFFICIENT FOR
C EACH SYMMETRY BLOCK

C-----------------------------------------------------------------------
CJDW  May/June 1996. Note on modifications of dimensioning.
C
C     Memory allocation :
C
C     It is assumed that there will never be more than (3*NATOM)**2
C     energy/gradient calculations. The memory allocation in the main
C     program for SETPTS and UPD_FD is based on this assumption.
C     Note that this bound is highly wasteful for GRADONLY calculations,
C     as the number of gradients will be much less. However, the code
C     currently reads record GRDPOINT even during ENERONLY jobs, hence
C     we need to make appropriate allowances.
C-----------------------------------------------------------------------

c INPUT
c integer NATOM
c integer NIRREP
c char*4  TYPE
c integer NDSCR

c OUTPUT
c integer NENER
c char*8  LABEL(NIRREP)
c integer ISYMIRR(3*NATOM)
c integer IPTTYPE(9*NATOM*NATOM)
c integer INVOP(3*NATOM)
c double  POINTS(27*NATOM*NATOM*NATOM)
c double  DSCR(NDSCR)

c RECORDS
c get 'NFDIRREP'
c get 'FDIRREP '
c get 'COORD   '
c get 'ATOMMASS'
c get TYPE//'SYMQ'
c get TYPE//'SYQT'
c get TYPE//'DEGN'
c get TYPE//'LABL'
c get 'INVPSMAT'
c get 'NUMVIBRT'
c put 'NUMPOINT'
c put 'FDCALCTP'
c put 'FDCOORDS'
c put 'NPTIRREP'
c put 'ENGPOINT'
c put 'GRDPOINT'
c put 'DIPPOINT'
c put 'POLPOINT'

      SUBROUTINE SETPTS(NATOM,NIRREP,TYPE,
     &                  NENER,LABEL,ISYMIRR,IPTTYPE,INVOP,
     &                  POINTS,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 TYPE
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ISYMIRR(3*NATOM),IPTTYPE(9*NATOM*NATOM),INVOP(3*NATOM)
      DIMENSION POINTS(27*NATOM*NATOM*NATOM)
      DIMENSION DSCR(NDSCR)

      INTEGER SKIPIR
      DIMENSION IJUNK(20),SKIPIR(20),idegen(100),nptirr(100)
      LOGICAL PRINTQ

      COMMON /FLAGS/  IFLAGS(100)
#include "control.com"
#include "machsp.com"
#include "parallel_aces.com"

      DATA TOL /1.D-8/

      PRINTQ=(IFLAGS(1).GE.10)
      NSIZE=3*NATOM

      if (gmtryopt.or.gradonly.or.nprocs.eq.1) then
         iRoot = 0
      else
c      o the root process is doing the reference geometry in a vib freq calc
c        with numerical gradients
         iRoot = 1
      end if

c   o check for restrictions on FD irreps
      CALL GETREC(-1,'JOBARC','NFDIRREP',1,ICOUNT)
      IF (ICOUNT.NE.0) THEN
         CALL GETREC(-1,'JOBARC','FDIRREP ',ICOUNT,IJUNK)
         WRITE(6,501)
501      FORMAT(T3,'@SETPTS: FCM evaluation limited ',
     &             'to the following symmetries :')
         WRITE(6,'((12I5))')(IJUNK(I),I=1,ICOUNT)
         DO I=1,20
            SKIPIR(I)=1
         END DO
         DO I=1,ICOUNT
            SKIPIR(IJUNK(I))=0
         END DO
      ELSE
         CALL IZERO(SKIPIR,20)
      END IF

      STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
      IF (PRINTQ) THEN
         WRITE(6,500)STPSIZ
500      FORMAT(T3,'Step size will be ',F8.5,' amu**(1/2) * bohr.')
      END IF
#ifdef _DEBUG_LVL0
      Print*, "The step size in symcor:", STPSIZ
#endif 
      lFREE = 1

      lSYMQ = lFREE
      lFREE = lFREE + NSIZE*NSIZE
      lCOOR = lFREE
      lFREE = lFREE + NSIZE
      lMASS = lFREE
      lFREE = lFREE + NATOM
      NDSCRLFT = NDSCR+1-lFREE
      IF (NDSCRLFT.LT.0) THEN
         print *, '@SETPTS: Insufficient memory.'
         print *, '         need ',-ndscrlft*ifltln,' more bytes'
         call aces_exit(1)
      END IF
      CALL GETREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE*IINTFP,
     &                                                DSCR(lSYMQ))
      CALL GETREC(20,'JOBARC','COORD   ',NSIZE*IINTFP,DSCR(lCOOR))
#ifdef _DEBUG_LVL0
      Print*, "The Cartesian coordinates read in septs:"
      Write(6,"(f10.4)"), (Dscr(I), I=lCOOR, NSIZE)
#endif      
      CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,DSCR(lMASS))

      CALL GETREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISYMIRR)
      CALL GETREC(20,'JOBARC',TYPE//'DEGN',NIRREP,IDEGEN)
      CALL GETREC(20,'JOBARC','INVPSMAT',3*NATOM,INVOP)
      IF (PRINTQ) THEN
      CALL GETREC(20,'JOBARC',TYPE//'LABL',NIRREP*IINTFP,LABEL)
      END IF

      CALL GETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      NLEFT=NMODE

      CALL ZERO(POINTS,27*NATOM*NATOM*NATOM)
      CALL IZERO(NPTIRR,NIRREP)

c   o set up vector of reciprocal square roots of atomic masses
      DO IOFF=0,NATOM-1
         X=SQRT(DSCR(lMASS+IOFF))
         IF (X.LT.TOL) THEN
            DSCR(lMASS+IOFF)=0.d0
         ELSE
            DSCR(lMASS+IOFF)=1.d0/X
         END IF
      END DO

c   o loop over irreducible representations and process first occurances
      NENER=0
      NGRAD=0
      NPTOTX=0
      NENERX=0
      NGRADX=0
      IOFF=0
      IPOS2=1
      IPOS3=1
      DO IRREP=1,NIRREP
         IFIRST=ISRCHEQ(NMODE,ISYMIRR,1,IRREP)
      IF (IFIRST.NE.NMODE+1) THEN

         ILAST=ISRCHNE(NLEFT,ISYMIRR(IFIRST),1,IRREP)
         NVIBSYM=ILAST-1
         NVIBUNQ=NVIBSYM/IDEGEN(IRREP)

         IF (PRINTQ) THEN
            WRITE(6,2000)LABEL(IRREP),IDEGEN(IRREP),NVIBUNQ
2000        FORMAT(T3,' Symmetry : ',A,' Degeneracy : ',I1,
     &             ' Unique symmetry coordinates : ',I3)
         END IF

         IPOS=NSIZE*(IFIRST-1)
#ifdef _DEBUG_LVL0
         Print*, "In setpts, energy only?", ENERONLY 
#endif
         IF (ENERONLY) THEN
            CALL DOENER(NATOM,NVIBUNQ,
     &                  DSCR(lSYMQ+IPOS),DSCR(lCOOR),STPSIZ,DSCR(lMASS),
     &                  POINTS(IPOS2),NPOINT,
     &                  INVOP(IPOS3),PRINTQ,DSCR(lFREE),NDSCRLFT)
            NENER=NENER+NPOINT
            IF (SKIPIR(IRREP).EQ.0) THEN
               do iR = 0, nProcs-1
                  call paces_batch_stat(iR,nProcs,iRoot,nPoint,iO,nE)
                  do iX = 1+iO, nE+iO
                     iPtType(iOff+iX) = 1+iR
                  end do
                  if (iR.eq.iRank) nEnerX=nEnerX+nE
               end do
               iRoot = mod(iRoot+nPoint,nProcs)
               nPTotX=nPTotX+nPoint
            ELSE
               DO IX=1,NPOINT
                  IPTTYPE(IOFF+IX)=0
               END DO
            END IF
         ELSE
#ifdef _DEBUG_LVL0
         Print*, "In setpts, grad only?, T"
#endif  

            CALL DOGRAD(NATOM,NVIBUNQ,
     &                  DSCR(lSYMQ+IPOS),DSCR(lCOOR),STPSIZ,DSCR(lMASS),
     &                  POINTS(IPOS2),NPOINT,
     &                  INVOP(IPOS3),PRINTQ,DSCR(lFREE),NDSCRLFT)
            NGRAD=NGRAD+NPOINT
            IF (SKIPIR(IRREP).EQ.0) THEN
               do iR = 0, nProcs-1
                  call paces_batch_stat(iR,nProcs,iRoot,nPoint,iO,nE)
                  do iX = 1+iO, nE+iO
                     iPtType(iOff+iX) = 1+iR
                  end do
                  if (iR.eq.iRank) nGradX=nGradX+nE
               end do
               iRoot = mod(iRoot+nPoint,nProcs)
               nPTotX=nPTotX+nPoint
            ELSE
               DO IX=1,NPOINT
                  IPTTYPE(IOFF+IX)=0
               END DO
            END IF
         END IF
         NPTIRR(IRREP)=NPTIRR(IRREP)+NPOINT
         IPOS2=IPOS2+NPOINT*NSIZE
         IPOS3=IPOS3+NVIBSYM
         IOFF=IOFF+NPOINT
         NLEFT=NLEFT-NVIBSYM

      END IF
      END DO

      NTOT=NGRAD+NENER
      NTOTX=NGRADX+NENERX
      WRITE(6,1000)NTOTX,NENERX,NGRADX
1000  FORMAT(T3,' Total number of calculations required      : ',I5,/,
     &       T3,' Number of single-point energy calculations : ',I5,/,
     &       T3,' Number of energy gradient     calculations : ',I5)
      IF (NTOT.GT.NSIZE*NSIZE) THEN
         WRITE(6,*) '@SETPTS: Too many points!'
         CALL ERREX
      END IF
      IF (NPTOTX.EQ.0) THEN
         WRITE(6,*) '@SETPTS: There are no vibrational modes!'
         CALL ERREX
      END IF
      IF (NTOTX.EQ.0) THEN
         if (nprocs.eq.1.or.irank.ne.0.or.gmtryopt.or.gradonly) then
            WRITE(6,*) '@SETPTS: There are no points to calculate!'
            CALL ERREX
         else
            CALL PUTREC(20,'JOBARC','LASTGEOM',1,1)
         end if
      END IF
C
C This record was added to do manual finite difference calculations
C with ACES III. The xjoda -procs n -rank m is run before execution
C is transfered to ACESIII. Ajith Perera, 10/2010.

      IF (.NOT. (nprocs.eq. 1 .and. irank.eq.0)) Call PUTREC(20, 
     &                               "JOBARC", "MANULFDS", 1, 1)


c   o write geometries to be used and calculation types to jobarc
      CALL PUTREC(20,'JOBARC','NUMPOINT',1,NTOT)
      CALL PUTREC(20,'JOBARC','FDCALCTP',NTOT,IPTTYPE)
      CALL PUTREC(20,'JOBARC','FDCOORDS',IINTFP*NSIZE*NTOT,POINTS)
      CALL PUTREC(20,'JOBARC','NPTIRREP',NIRREP,NPTIRR)

c   o write out zero vectors for calculation archiving
      CALL ZERO(POINTS,NTOT*max(NSIZE,9))
      CALL PUTREC(20,'JOBARC','ENGPOINT',IINTFP*      NTOT,POINTS)
      CALL PUTREC(20,'JOBARC','GRDPOINT',IINTFP*NSIZE*NTOT,POINTS)
      CALL PUTREC(20,'JOBARC','DIPPOINT',IINTFP*    3*NTOT,POINTS)
      CALL PUTREC(20,'JOBARC','POLPOINT',IINTFP*    9*NTOT,POINTS)
      

      RETURN
      END

