      SUBROUTINE GENTRANMAT(SCR, MAXCOR)
C
C MAXPRIM  - Maximum number of primitives in a shell (30).
C MAXFNC  - as above (apparently). Number of contracted functions in
C           a shell. Sometimes MAXPRIM is used instead.
C           [MAXPRIM and MAXFNC must keep same value the way things
C                are coded. It is desirable in any case to have them the
C                same so we can do large uncontracted calculations.
C                However, the allocation for READIN and dimensioning in
C                READIN are a bit sloppy, and should be improved.]
C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers (400).
C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers (800).
C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers (100).
C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions) (500).
C
C     /INDX/  :
C
C     KMAX          - Number of shells.
C     KHKT(7)       - KHKT(I) = I*(I+1)/2
C     NHKT(MXTNSH)  - L+1 value for each shell (1 for s, 2 for p, etc).
C     NUCO(MXTNSH)  - Number of primitives in each shell.
C     NRCO(MXTNSH)  - Number of contracted functions in each shell.
C
C     /DAT/   :
C
C     ALPHA(MXTNPR)   - Exponents of symmetry inequivalent centers.
C     CONT(MXTNCC)    - Coefficients of symmetry inequivalent centers.
C     CENT(3,MXTNSH ) - Coordinates of symmetry inequivalent shells.
C     CORD(100,3)     - Coordinates of symmetry inequivalent centers.
C     CHARGE(100)     - Nuclear charges of symmetry inequivalent centers.
C     FMULT(8)        - 
C     TLA             - A cutoff 10**-INTGRL_TOL.
C     TLC             - 0.1*TLA ?
C
C     /VMTASK/ :
C     
C     ITASK           -  0  NORMAL CALCULATION (DEFAULT)
C                     -  1  READIN ONLY
C                     -  2  READIN+ONEL+ONELH ONLY
C                     - >3  NORMAL CALCULATION.
C
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Common block declarations (or at least some).
C-----------------------------------------------------------------------
      INTEGER NEWIND,MSTOLD
      DOUBLE PRECISION PC
      INTEGER DSTRT,NTAP,LU2,NRSS,NUCZ,ITAG,MAXLOP,MAXLOT,KMAX,NMAX,
     &        KHKT,MULT,ISYTYP,ITYPE,AND,OR,EOR,NPARSU,NPAR,MULNUC,
     &        NHKT,MUL,NUCO,NRCO,JSTRT,NSTRT,MST,JRS
      DOUBLE PRECISION ALPHA,CONT,CENT,CORD,CHARGE,FMULT,TLA,TLC
      DIMENSION SCR(MAXCOR)
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
      parameter (mxp2=maxprim*maxprim)
      parameter (khm=(nht*(nht+1))/2,nhl=mxp2*khm*khm)
      parameter (kwd=(nht+1)*(nht+2)*(nht+3)/6,nh4=4*nht-3)
C
      COMMON /TST/ TIM(40),IFREQ(40)
      COMMON /MMMM/ M2(12)
      COMMON /REP/ NEWIND(MXCBF) ,MSTOLD(8)
      COMMON /FLAGS/ IFLAGS(100)
C
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     & MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     & AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(100),
     & NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     & NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
C
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &              CORD(100,3),CHARGE(100),FMULT(8),TLA, TLC
C
      COMMON /SYMIND/ IBFS(MXCBF)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /MEMINF/ IJUNK(2),ITOTMEM,IHWM1,IHWM2
C-----------------------------------------------------------------------
C
      IPRINT = IFLAGS(1)
      INTES=0
      INTE53=0
      INTE32=0
      INTEF=4
C
      CALL ZERO(TIM,40)
      CALL IZERO(IFREQ,40)
      CALL  INITP
C
C ALLOCATE CORE FOR READIN ROUTINE
C
      I000 = 1
      I010 = I000 + MAXATM
      I020 = I010 + MXCBF*224
      I030 = I020 + MXTNCC
      I040 = I030 + MXTNCC
      I050 = I040 + MAXFNC*MAXFNC
      I060 = I050 + MAXFNC
      I070 = I060 + MAXFNC
      I080 = I070 + MAXFNC
      I090 = I080 + MXSHEL*28*15*8
      I100 = I090 + 50
      I110 = I100 + KWD*3
      I120 = I110 + NHT*MAXFNC
      I130 = I120 + NHT*MAXFNC
      I140 = I130 + 8
      I150 = I140 + 100*3
      I160 = I150 + 8*3
      I170 = I160 + MXTNCC
      I180 = I170 + NHT
      I190 = I180 + 8
      I200 = I190 + 18*8
      I210 = I200 + MXCBF
      I220 = I210 + MAXFNC*2
      I230 = I220 + MXCBF*224
      I240 = I230 + MXTNPR
      I250 = I240 + 100
      I260 = I250 + 100
      I270 = I260 + 3
      I280 = I270 + MXCBF
      I290 = I280 + MXCBF*224
      I300 = I290 + MXCBF*224
      ITOP = I300 + MXCBF*224
C
      IF(ITOP .GE. MAXCOR)THEN
         CALL INSMEM("GENTRANMAT", ITOP, MAXCOR)
      ENDIF
C
      CALL BLTAOTOSO(Scr(I000),Scr(I010),Scr(I020),Scr(I030),Scr(I040),
     &               Scr(I050),Scr(I060),Scr(I070),Scr(I080),Scr(I090),
     &               Scr(I100),Scr(I110),Scr(I120),Scr(I130),Scr(I140),
     &               Scr(I150),Scr(I160),Scr(I170),Scr(I180),Scr(I190),
     &               Scr(I200),Scr(I210),Scr(I220),Scr(I230),Scr(I240),
     &               Scr(I250),Scr(I260),Scr(I270),MAXFNC,
     &               NHT,KWD,NH4,Scr(I280),Scr(I290),Scr(I300))
C
      RETURN
      END
