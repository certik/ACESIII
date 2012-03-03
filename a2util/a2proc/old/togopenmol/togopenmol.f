      PROGRAM TOGOPENMOL
c     
C     This program generates a cube file from an ACESII
C     calculation.
C     
C     Written by Tommy Miller, Quantum Theory Project,
C     University of Florida, 1998
C     
C     A portion of the source code and a number
C     of the subroutines were extracted from the HFDFT code
C     written by Nevin Oliphant, Quantum Theory Project,
C     University of Florida, 1993. Also, I would like to
C     thank Ken Wilson and Drs. John Watts and Ajith Perera for their
C     valuable time and assistance.
c     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
C     Characters for determining symmetry elements
      CHARACTER*1 SYMEL(3,3)
      CHARACTER*8 PGRP
      CHARACTER*12 STRIAORBO
      CHARACTER*12 STRIAORBD
      CHARACTER*12 STRICORR
      CHARACTER*12 STRIDORO
      CHARACTER*12 STRINORBS
      CHARACTER*12 STRINUMCUBE
      integer err
C     
      LOGICAL::YESNO
C     
      DIMENSION NOCC(16),ISHL(6),IORBVEC(9)
C     
      COMMON//ICORE(1)
      COMMON /ISTART/ I0
      COMMON /IOPOS/ ICRSIZ,ICHSZ,IOFF(2),LENREC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(10)
      COMMON /FLAGS2/ IFLAGS2(50)
      COMMON /ISYMINF/ NIRREP,NSOIRP(8)
      COMMON /CNST/ALPHA,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15
      COMMON /PAR/ PI,PIX4
      COMMON /IPAR/ LUOUT
C     
CXXXXXXXXXXXXXXXXXXXXX
C     READ IN THE OPTIONS FROM THE ZMAT FILE
CXXXXXXXXXXXXXXXXXXXXX
C     
C     Open the ZMAT file and search for the cube_gen namelist
      CALL NL_INIT('cube_gen',err,.true.)
      if(err.eq.1)then
         write(*,*)' *cube_gen namelist not found in ZMAT file'
         call errex
      end if
C     
C     Determine the number of cube files to be created.
      CALL NL_STR('NUM_CUBE','ONE',STRINUMCUBE)
C     
C     Interpret STRINUMCUBE:
      IF (STRINUMCUBE.EQ.'ONE') THEN
         NUMCUBE=1
         IDDIFF=0
      ELSE IF (STRINUMCUBE.EQ.'TWO') THEN
         NUMCUBE=2
         IDDIFF=0
      ELSE IF (STRINUMCUBE.EQ.'DEN_DIFF') THEN
         NUMCUBE=2
         IDDIFF=1
      END IF
C     
C     Loop over the cube files.
C     
      DO 911 ICUBE=1,NUMCUBE
C     
         IF (ICUBE.EQ.2) CALL NL_INIT('cube2_gen',err,.true.)
C     
C     Read in the upper and lower boundaries of the cube
         CALL NL_REAL('TVALX', 3.0D+00,TVALX)
         CALL NL_REAL('BVALX',-3.0D+00,BVALX)
         CALL NL_REAL('TVALY', 3.0D+00,TVALY)
         CALL NL_REAL('BVALY',-3.0D+00,BVALY)
         CALL NL_REAL('TVALZ', 3.0D+00,TVALZ)
         CALL NL_REAL('BVALZ',-3.0D+00,BVALZ)
C     
C     Read in the number of interval along each axis
         CALL NL_INT('NPTSX',27,NPTSX)
         CALL NL_INT('NPTSY',27,NPTSY)
         CALL NL_INT('NPTSZ',27,NPTSZ)
     
C     Determine if this is a OFT calculation 
c     (this keyword is not for the masses yet)
         CALL NL_INT('DFT',0,NDFT)
C     Read in the various options regarding what is to be calculated
         CALL NL_STR('CUBE',      'DENSITY',STRIDORO)
         CALL NL_STR('DEN_METHOD','SCF',    STRICORR)
         CALL NL_STR('DENSITY',   'TOTAL',  STRIAORBD)
         CALL NL_STR('ORBITAL',  'ALPHA',  STRIAORBO)
         CALL NL_STR('ORB_METHOD','SCF',    STRINORBS)
C     
C     Translate string options
C     
C     If IDORO=0, just density is displayed
C             =1, just orbitals are displayed
C             =2, both are displayed
         IF (STRIDORO.EQ.'DENSITY') IDORO=0
         IF (STRIDORO.EQ.'ORBITALS') IDORO=1
         IF (STRIDORO.EQ.'BOTH') IDORO=2
C     
C     If ICORR=0, just SCF density is displayed
C             =1, SCF+relaxation density is displayed (only for
C     correlated methods)
         IF (STRICORR.EQ.'SCF') ICORR=0
         IF (STRICORR.EQ.'CORRELATED') ICORR=1
C     
C     If INORBS=0, SCF orbitals will be displayed
C              =1, Correlated Natural Orbitals will be displayed
         IF (STRINORBS.EQ.'SCF') INORBS=0
         IF (STRINORBS.EQ.'CORRELATED') INORBS=1
C     
C     IF IAORBD=0, just alpha density is displayed
C              =1, just beta density is displayed
C              =2, total density is displayed
         IF (STRIAORBD.EQ.'ALPHA') IAORBD=0
         IF (STRIAORBD.EQ.'BETA')  IAORBD=1
         IF (STRIAORBD.EQ.'TOTAL') IAORBD=2
C     
C     IF IAORBO=0, just alpha orbitals are displayed
C              =1, just beta orbitals are displayed
C              =2, both are displayed
C              =3, certain selected orbitals are displayed
         IF (STRIAORBO.EQ.'ALPHA') IAORBO=0
        IF (STRIAORBO.EQ.'BETA')  IAORBO=1
        IF (STRIAORBO.EQ.'BOTH')  IAORBO=2
        IF (STRIAORBO.EQ.'SELECT')  THEN
           IAORBO=3
           CALL NL_INT('ORBITAL1',0,IORBVEC(1))
           CALL NL_INT('ORBITAL2',0,IORBVEC(2))
           CALL NL_INT('ORBITAL3',0,IORBVEC(3))
           CALL NL_INT('ORBITAL4',0,IORBVEC(4))
           CALL NL_INT('ORBITAL5',0,IORBVEC(5))
           CALL NL_INT('ORBITAL6',0,IORBVEC(6))
           CALL NL_INT('ORBITAL7',0,IORBVEC(7))
           CALL NL_INT('ORBITAL8',0,IORBVEC(8))
           CALL NL_INT('ORBITAL9',0,IORBVEC(9))
        END IF
        call nl_term
C     
C     When Natural Orbitals are requested, the correlated is exclusively calculated
        IF (INORBS.EQ.1) THEN
           ICORR=1
           WRITE(*,*) "The correlated density is"
           WRITE(*,*) "exclusively calculated with NOs."
        END IF
C     
C     WRITE(*,*) "IDORO=",IDORO
C     WRITE(*,*) "ICORR=",ICORR
C     WRITE(*,*) "IAORBD=",IAORBD
C     WRITE(*,*) "IAORBO=",IAORBO
C     WRITE(*,*) "INORBS=",INORBS
C     
CXXXXXXXXXXXXXXXXXXX
C     
        PI=DATAN(1.0D+00)*4.D+00
        PIX4=DATAN(1.0D+0)*16.D+00
        LUOUT=6
C     
C     INITIALIZE THE PROGRAM BY CALLING CRAPSI, GET IUHF FLAG
C     
        IF (ICUBE.EQ.1) THEN
C     
           CALL CRAPSI(ICORE(1),IUHF,0)
C     
C     MAXIMUM CORE AVAILABLE FOR THIS CALCULATION
C     
           CALL GETREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
           MAXCOR=ICRSIZ
C     
C     Output control
           IPRINT=IFLAGS(1)
C     Print only results
           IF(IPRINT.EQ.1) IPRINT=0
C     Print intermediate information
           IF(IPRINT.GT.1) IPRINT=1
C     
           IF (IPRINT.EQ.1) THEN
              IF (IUHF.EQ.1) THEN
                 WRITE(LUOUT,1010)
              ELSE
                 WRITE(LUOUT,1020)
              ENDIF
           ENDIF
 1010      FORMAT(/' THE CALCULATION IS UHF')
 1020      FORMAT(/' THE CALCULATION IS RHF')
C     
C     
C     Get the total number of centers
           CALL GETREC(20,'JOBARC','NATOMS  ',1,NATOMS)
C     
C     
C     Get the number of unique atoms, (orbits), in the full point group
           CALL GETREC(20,'JOBARC','FULLNORB',1,IUATMS) 
C     
C Get the number of unique atoms, (orbits), in the computational
C     point group
           CALL GETREC(20,'JOBARC','COMPNORB',1,IUCATMS)
        END IF
C     
C     The atomic number of each center
        INUC=I0
C     The coordinates of each center
        ICOORD=INUC+NATOMS + mod(natoms,2)
C     The number of atoms in each orbit in the computational point group
        IPOPC=ICOORD+3*NATOMS*IINTFP
C     The atoms sorted by computational point group orbits
        IMEMBC=IPOPC+IUCATMS + mod(IUCATMS,2)
C     
C     Starting point of remaining memory
        ISTART=IMEMBC+NATOMS + mod(natoms,2)
C     
C     Get the number of atoms in a computational point group orbit
        CALL GETREC(20,'JOBARC','COMPPOPV',IUCATMS,ICORE(IPOPC))
C     Get the sorted list of atoms for the computational point group
        CALL GETREC(20,'JOBARC','COMPMEMB',NATOMS,ICORE(IMEMBC))
C     
C     Get the atomic numbers and coordinates
        CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS,ICORE(INUC))
        CALL GETREC(20,'JOBARC','COORD',3*IINTFP*NATOMS,ICORE(ICOORD))
C     
C     Get the full point group
        CALL GETREC(20,'JOBARC','FULLPTGP',IINTFP,PGRP)
C     
        WRITE(LUOUT,3010) PGRP
 3010   FORMAT(/'The Full Point Group is ',A4)
C     
C     Get the computational point group
        CALL GETREC(20,'JOBARC','COMPPTGP',IINTFP,PGRP)
C     
        WRITE(LUOUT,3020) PGRP
 3020   FORMAT(/' The Computational Point Group is ',A4)
C     
        IF (IPRINT.EQ.1) THEN
           WRITE(LUOUT,1102) IUCATMS
        ENDIF
 1102   FORMAT(' THERE ARE',I4,
     $       ' UNIQUE ATOMS IN THE COMPUTATIONAL POINT GROUP')
C     
C     Get number of basis functions
        CALL GETREC(20,'JOBARC','NAOBASFN',1,NBAS)
C     
        IF (IPRINT.EQ.1)THEN
           WRITE(LUOUT,1030) NBAS
        ENDIF
 1030   FORMAT(/' THE NUMBER OF ATOMIC BASIS FUNCTIONS IS',I5)
C     
C     
C     Get number of irreducible representations.
        CALL GETREC(20,'JOBARC','COMPNIRR',1,NIRREP)
C     Get number of symmetry orbitals in each irrep
        CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,NSOIRP)
C     
C     Get alpha and beta occupation
        CALL GETREC(20,'JOBARC','OCCUPYA',NIRREP,NOCC(1))
        NOCCA=0
        DO 10 I=1,NIRREP
           NOCCA=NOCCA+NOCC(I)
 10     CONTINUE
        IF (IUHF.EQ.1) THEN
           CALL GETREC(20,'JOBARC','OCCUPYB',NIRREP,NOCC(9))
           NOCCB=0
           DO 20 I=9,NIRREP+8
              NOCCB=NOCCB+NOCC(I)
 20        CONTINUE
        ELSE
           NOCCB=NOCCA
        ENDIF
C     
        IF (IPRINT.EQ.1)THEN
           WRITE(LUOUT,1040) NOCCA
           WRITE(LUOUT,1050) NOCCB
        ENDIF
 1040   FORMAT(/' THE ALPHA OCCUPATION IS',I4)
 1050   FORMAT(' THE BETA OCCUPATION IS',I4)
C     
C     First get the atomic numbers of the atoms, the number of
C     primitive functions for each unique atom, the coordinates for
C     each unique atom, the total number of primitive functions for
C     the molecule, the total number of primitive coefficients for the
C     molecule, the number of different angular momenta for each unique
C     atom (NUMOM(IUCATMS), the largest number of angular momenta for a
C     single atom (NTANGM).
C     
C     The number of primitive functions for each unique atom
        NUFCT=ISTART
C     The number of angular momentum shells for each unique atom
        NUMOM=NUFCT+IUCATMS + mod(IUCATMS,2)
C     The number of primitive functions for each atom
        NFCT=NUMOM+IUCATMS + mod(IUCATMS,2)
C     The number of angular momentum shells for each atom
        NANGMOM=NFCT+NATOMS + mod(NATOMS,2)
C     Array of atomic names
        IATMNAM=NANGMOM+NATOMS + mod(NATOMS,2)
C     Number of AOs for each unique atom
        NAOUATM=IATMNAM+IUCATMS*IINTFP
C     Number of AOs for each atom
        NAOATM=NAOUATM+IUCATMS + mod(IUCATMS,2) 
C     
      ISTART=NAOATM+NATOMS + mod(NATOMS,2)
C     
      IF(ISTART.GE.MAXCOR) CALL INSMEM('HYBRD1',ISTART,MAXCOR)
C     
      CALL BASIS(IUCATMS,NATOMS,ITFCT,LNP1,LNPO,NTANGM,ICORE(IMEMBC),
     $     ICORE(INUC),ICORE(NFCT),ICORE(NUFCT),ICORE(NANGMOM),
     $     ICORE(NUMOM),ICORE(IATMNAM),ICORE(ICOORD),ICORE(IPOPC),
     $     ICORE(NAOATM),ICORE(NAOUATM),IPRINT,ISHL)
C     
      write(*,*)" number of primitives for each atom"
      call printivec(icore(nfct),natoms)
      write(*,*)" largest number of primitives in a single shell: "
     $     ,lnp1
      write(*,*)" largest number of primitive orbitals in a ",
     $     "single shell: ",lnpo 
      write(*,*)" total number of primitive functions: ",itfct
      write(*,*)" number of AOs for each atoms:"
      call printivec(icore(naoatm),natoms)
c
C     Now fill the ao to primitive transformation, alpha and angular
C     momentum matrices.
C     
C     Get number of symmetry adapted basis functions
C     is different from nbas if spherical=on
      IF (IFLAGS(2).NE.0) THEN
         CALL GETREC(20,'JOBARC','NUMDROPA',1,NDROP)
         WRITE(*,*) "Ndrops=",NDROP
         IF (NDROP.NE.0) THEN
            WRITE(*,*) "NO FROZEN CORE!!!!!!."
         END IF
      END IF
      CALL GETREC(20,'JOBARC','NBASTOT ',1,NBASP)
      WRITE(*,*) "# orbitals =",NBASP
C     
      NVRTA=NBASP-NOCCA
      NVRTB=NBASP-NOCCB
C     
C     Note that if IAORBO=3, nbaspeff will redefined later.
      IF (IAORBO.EQ.2) THEN
         NBASPEFF=NBASP*2
      ELSE
         NBASPEFF=NBASP
      END IF
C     
C     The number of primitive functions for each angular momentum
C     shell in each atom
      NMOMFCT=ISTART
C     The number of atomic orbitals for each angular momentum
C     shell in each atom
      NMOMAO=NMOMFCT+NATOMS*NTANGM + mod(NATOMS*NTANGM,2)
C     The alpha matrix, (exponential factor).
         IALPHA=NMOMAO+NATOMS*NTANGM + mod(NATOMS*NTANGM,2)
C     
C     C1 MO coefficients
         ICOEFFA=IALPHA+ITFCT*IINTFP
         ICOEFFB=ICOEFFA+NBASP*NBAS*IINTFP
C     Natural Orbitals in the MO basis
         IANOCOEFF=ICOEFFB+NBASP*NBAS*IINTFP
         IBNOCOEFF=IANOCOEFF+NBASP*NBASP*IINTFP
C     C1 density matrices
         IDENSA=IBNOCOEFF+NBASP*NBASP*IINTFP
         IDENSB=IDENSA+NBASP*NBASP*IINTFP
C     Primitive coefficient matrix
         IPCOEFFA=IDENSB+NBASP*NBASP*IINTFP
         IPCOEFFB=IPCOEFFA+NBASP*ITFCT*IINTFP
C     AO-basis to primitive
         IPCOEFF=IPCOEFFB+NBASP*ITFCT*IINTFP
         IPCOEFF2=IPCOEFF+ITFCT*IINTFP*NBAS
C     Coordinates of the Array of points in cube
         IPCOORD=IPCOEFF2+ITFCT*IINTFP*NBAS
C     C1 primitive density matrices
         IAODENSA=IPCOORD+(NPTSZ+1)*3*IINTFP
         IAODENSB=IAODENSA+NBAS*NBAS*IINTFP
C     Total value of a,b,and tot density at each point
         IDENSACUB=IAODENSB+NBAS*NBAS*IINTFP
         IDENSBCUB=IDENSACUB+(NPTSZ+1)*IINTFP
         IDENSTCUB=IDENSBCUB+(NPTSZ+1)*IINTFP
C     Total value of each mo at each point
         IMOACUB=IDENSTCUB+(NPTSZ+1)*IINTFP
         IMOBCUB=IMOACUB+(NPTSZ+1)*IINTFP*NBASP
         IMOOUT=IMOBCUB+(NPTSZ+1)*IINTFP*NBASP
C     SALC to primitive tranformation matrix
         IPCOEFSO=IMOOUT+(NPTSZ+1)*IINTFP*NBASPEFF
C     Scratch space for reading MOL file and renormalization
         ISCR1=IPCOEFSO+ITFCT*IINTFP*NBASP
         ISCR2=ISCR1+LNP1*IINTFP
         ISTOP=ISCR2+LNPO*IINTFP
         ISTART=ISCR1
C     
         IF(ISTOP.GE.MAXCOR) CALL INSMEM('HYBRD2',ISTOP,MAXCOR)
C     
         CALL PRIM(NATOMS,IUCATMS,ITFCT,NBAS,LNP1,LNPO,NTANGM,
     $        ICORE(IPOPC),ICORE(NFCT),
     $        ICORE(NANGMOM),ICORE(NMOMFCT),ICORE(NMOMAO),ICORE(IMEMBC),
     $        ICORE(IALPHA),ICORE(IPCOEFF),ICORE(NAOATM),ICORE(ISCR1),
     $        ICORE(ISCR2),ISHL,MAXANG)
C
         write(*,*)" gaussian exponents"
         call printdvec(icore(iscr1),lnp1)
         write(*,*)" contraction coefficents"
         call printdvec(icore(iscr2),lnpo)
c
         IF ((IDORO.NE.1).OR.(INORBS.EQ.1)) THEN
            ICOEF=ISTART
            IDENRELA=ICOEF+NBASP*NBASP*IINTFP
            IDENRELB=IDENRELA+NBASP*NBASP*IINTFP
            IDENSCF=IDENRELB+NBASP*NBASP*IINTFP
            ICOEF2=IDENSCF+NBASP*NBASP*IINTFP
            ITEMPMAT=ICOEF2+NBASP*NBASP*IINTFP
            IREORD=ITEMPMAT+NBASP*NBASP*IINTFP
            INOCOEFFTEMP=IREORD+NBASP+MOD(NBASP,2)
            INOOCC=INOCOEFFTEMP+NBASP*NBASP*IINTFP
            IWORK=INOOCC+NBASP*NBASP*IINTFP
            IUMAT=IWORK+(3*NBASP-1)*IINTFP
            ISTOP=IUMAT+NBAS*NBASP*IINTFP
C     
            IF(ISTOP.GE.MAXCOR) CALL INSMEM('HYBRD3',ISTART,MAXCOR)
C     
C     Get scf density matrix
            IFLG=0
            CALL MKDEN(ICORE(ICOEF),ICORE(IDENSA),NBASP,NOCC,IFLG,ICORR,
     $           ICORE(IDENRELA),ICORE(IDENSCF),ICORE(ICOEF2),
     $           ICORE(ITEMPMAT),IUHF,NDFT)
            IF (IUHF.EQ.1)THEN
               IFLG=1
               CALL MKDEN(ICORE(ICOEF),ICORE(IDENSB),NBASP,NOCC,IFLG,
     $              ICORR,
     $              ICORE(IDENRELB),ICORE(IDENSCF),ICORE(ICOEF2),
     $              ICORE(ITEMPMAT),IUHF,NDFT)
            ELSE
               IFLG=1
               CALL MKDEN(ICORE(ICOEF),ICORE(IDENSB),NBASP,NOCC,IFLG,
     $              ICORR,
     $              ICORE(IDENRELB),ICORE(IDENSCF),ICORE(ICOEF2),
     $              ICORE(ITEMPMAT),IUHF,NDFT)
            ENDIF
C     
            IF (IDORO.NE.1) THEN
C     Form the SALC to primitive transformation matrix
C     
C     Get the SALC to AO transformation matrix
               CALL GETREC(20,'JOBARC','CMP2ZMAT',NBAS*NBASP*IINTFP,
     $              ICORE(IUMAT))
C     Left multiply by the prim to AO matrix to get the sale to prim. mat.
               CALL XGEMM('N','N',ITFCT,NBASP,NBAS,1.0D+00,
     $              ICORE(IPCOEFF2),
     $              ITFCT,ICORE(IUMAT),NBAS,0.0D+00,ICORE(IPCOEFSO),
     $              ITFCT)
C     
               WRITE(*,*)" AO to prim #1"
               CALL PRINTMAT(ICORE(IPCOEFF),ITFCT,NBAS)
               WRITE(*,*)" AO to prim #2"
               CALL PRINTMAT(ICORE(IPCOEFF2),ITFCT,NBAS)
               WRITE(*,*)" salc to ao"
               CALL PRINTMAT(ICORE(IUMAT),NBAS,NBASP)
               WRITE(*,*)" Rsalc to prim transformation"
               CALL PRINTMAT(ICORE(IPCOEFSO),ITFCT,NBASP)
            END IF
C     
            IF (INORBS.EQ.1) THEN
               CALL GETREC(20,'JOBARC','REORDERA',NBASP,ICORE(IREORD))
               WRITE(*,*) "Reorder vector"
               CALL PRINTIVEC(ICORE(IREORD),NBASP)
C     
C     Diagonalize the alpha and beta correlated density
C     to obtain the NOs in an MO basis
C     Important Output: IANOCOEFF and IBNOCOEFF
C     
               CALL NORBMAKE(ICORE(IDENRELA),ICORE(IDENRELB),
     $              ICORE(INOCOEFFTEMP),ICORE(IANOCOEFF),
     $              ICORE(IBNOCOEFF),
     $              ICORE(IREORD),NBASP,ICORE(INOOCC),
     $              ICORE(IWORK))
            END IF
         END IF
C     
C     Alpha
         IFLG=0
         IFLG2=0
C     Scratch arrays for sorting MO vectors
         ICNTR=ISTART
         ICOEFF=ICNTR+NBAS + mod(nbas,2)
         ICOEFFTEMP=ICOEFF+NBASP*NBAS*IINTFP
         ISTOP=ICOEFFTEMP+NBASP*NBAS*IINTFP
         ISTART=ICNTR+NBAS + mod(nbas,2)
C     
         IF(ISTOP.GE.MAXCOR) CALL INSMEM('HYBRD5',ISTOP,MAXCOR)
C     
         IF (IDORO.NE.0) THEN
C Get MO vectors in correct order in C1 symmetry
C     
            CALL GTAO(IFLG,IFLG2,NATOMS,NBAS,NBASP,NTANGM,ICORE(ICNTR),
     $           ICORE(ICOEFF),ICORE(ICOEFFA),ICORE(NANGMOM),
     $           ICORE(NMOMAO),ICORE(NAOATM),IPRINT,ICORE(IANOCOEFF),
     $           INORBS,ICORE(ICOEFFTEMP))
C     
C     Form the C1 primitive coefficient vector
            CALL XGEMM('N','N',ITFCT,NBASP,NBAS,1.0D+00,
     $           ICORE(IPCOEFF),
     $           ITFCT,ICORE(ICOEFFA),NBAS,0.0D+00,
     $           ICORE(IPCOEFFA),ITFCT)
C     
         END IF
C     
C     Scratch arrays for primitive values at a point in the cube
         IPVAL=ISTART
         IPVAL2=IPVAL+ITFCT*NBASP*IINTFP
         IPCRDX=IPVAL2+ITFCT*NBASP*IINTFP
         IPCRDY=IPCRDX+(NPTSX+1)*IINTFP
         IPCRDZ=IPCRDY+(NPTSY+1)*IINTFP
         IIPCNT=IPCRDZ+(NPTSZ+1)*IINTFP
         IAPDX=IIPCNT+NBAS+MOD(NBAS,2)
         IAPDY=IAPDX+NATOMS*IINTFP
         IAPDZ=IAPDY+NATOMS*IINTFP
         IAPDR=IAPDZ+NATOMS*IINTFP
         IDENSNMAT=IAPDR+NATOMS*IINTFP
         IFINMO=IDENSNMAT+NBASP*NBASP*IINTFP
         IFINSALC=IFINMO+NBASP*IINTFP
         INUC2=IFINSALC+NBASP*IINTFP
         ICOORD2=INUC2+NATOMS+mod(NATOMS,2)
         IIORBVECA=ICOORD2+3*NATOMS*IINTFP
         IIORBVECB=IIORBVECA+10
         IORDERMATA=IIORBVECB+10
         IORDERMATB=IORDERMATA+9*IINTFP*ITFCT
         ISTOP=IORDERMATB+9*IINTFP*ITFCT
C     
         IF(ISTOP.GE.MAXCOR) CALL INSMEM('HYBRD6',ISTOP,MAXCOR)
C     
C     
C     If needed, reduce the primitive coeff. matrices to only contain
C     the desired orbitals.
         IF (IAORBO.EQ.3) THEN
            CALL ORBOUT(NBASP,NORBOUTA,NORBOUTB,ICORE(IPCOEFFA),
     $           ICORE(IPCOEFFB),ICORE(IIORBVECA),ICORE(IIORBVECB),
     $           ICORE(IORDERMATA),ICORE(IORDERMATB),ITFCT,NORBOUT,
     $           IORBVEC)
            NBASPEFF=NORBOUT
         ELSE
            NORBOUT=NBASPEFF
            NORBOUTA=NBASP
            NOBROUTB=NBASP
         END IF
C
C     WRITE(*,*) "# of ang. mom. shells for each atom"
C     CALL PRINTIVEC(ICORE(NANGMOM),NATOMS)
C     WRITE(*,*) "# of prim. funcs. for each ang.
C     $ mom. shell in each atom"
C     CALL PRINTIVEC(ICORE(NMOMFCT),(NATOMS*NTANGM))
C     WRITE(*,*) "# of AOs for each atom"
C     CALL PRINTIVEC(ICORE(NAOATM),NATOMS)
C     WRITE(*,*) "Alpha coeffs."
C     CALL PRINTDVEC(ICORE(IALPHA),ITFCT)
C     WRITE(*,*) "gaussian coeffs."
C     CALL PRINTMAT(ICORE(IPCOEFFA),ITFCT,NORBOUTA)
C     
C     WRITE(*,*) "atoms"
C     CALL PRINTMAT(ICORE(ICOORD),3,NATOMS)
C     
C     Heading for output
C     
         NATMS=NATOMS
         CALL REMOVE(NATOMS,NATMS,ICORE(INUC),ICORE(ICOORD),
     $        ICORE(ICOORD2),ICORE(INUC2))
C     
         IF (ICUBE.EQ.1) THEN
            IF (IDORO.NE.1) OPEN(UNIT=98,FILE='den.cube')
            IF (IDORO.NE.0) OPEN(UNIT=99,FILE='orb.cube')
            IDENOUT=98
            IORBOUT=99
         ELSE
            IF (IDORO.NE.1) OPEN(UNIT=96,FILE='den2.cube')
            IF (IDORO.NE.0) OPEN(UNIT=97,FILE='orb2.cube')
            IDENOUT=96
            IORBOUT=97
         END IF
C     
         CALL HEADER(NATMS,ICORE(INUC2),ICORE(ICOORD2),NPTSX,
     $        TVALX,BVALX,NPTSY,TVALY,BVALY,NPTSZ,TVALZ,BVALZ,
     $        NBASPEFF,IDORO,IAORBD,ICORR,INORBS,ICORE(IIORBVECA),
     $        ICORE(IIORBVECB),NORBOUTA,NORBOUTB,IAORBO,NBASP,
     $        IORBOUT,IDENOUT)
C     
         DO 100 IBIGLOOP=1,(NPTSX+1)
            DO 101 JBIGLOOP=1,(NPTSY+1)
C     Determine the grid points
               CALL GRID(NPTSX,TVALX,BVALX,NPTSY,TVALY,BVALY,NPTSZ,
     $              TVALZ,BVALZ,ICORE(IPCOORD),ICORE(IPCRDX),
     $              ICORE(IPCRDY),ICORE(IPCRDZ),IBIGLOOP,JBIGLOOP)
C     WRITE(*,*) "points in row ",((IBIGLOOP-l)*(NPTSY+l)+JBIGLOOP)
C     CALL PRINTMAT(ICORE(IPCooRD),3,(NPTSZ+l))
C     
C     Determine the desired quantities at the grid points for alpha
               CALL GRIDVAL(ICORE(IPCOORD),NPTSZ,ICORE(IPVAL),ITFCT,
     $              NBASP,
     $              ICORE(ICOORD),NATOMS,ICORE(IALPHA),
     $              ICORE(IPCOEFFA),
     $              ICORE(NMOMFCT),MAXANG,NTANGM,NBAS,
     $              ICORE(IDENSA),ICORE(IDENSACUB),0,ICORE(IMOACUB),
     $              ICORE(IPVAL2),ICORE(IIPCNT),ICORE(IAPDX),
     $              ICORE(IAPDY),ICORE(IAPDZ),ICORE(IAPDR),
     $              ICORE(IDENSNMAT),
     $              ICORE(IPCOEFSO),ICORE(IFINMO),ICORE(IFINSALC),
     $              IBIGLOOP,
     $              JBIGLOOP,IDORO,NORBOUTA)
C     Determine the desired quantities at the grid points for beta
               CALL GRIDVAL(ICORE(IPCOORD),NPTSZ,ICORE(IPVAL),ITFCT,
     $              NBASP,
     $              ICORE(ICOORD),NATOMS,ICORE(IALPHA),
     $              ICORE(IPCOEFFB),
     $              ICORE(NMOMFCT),MAXANG,NTANGM,NBAS,
     $              ICORE(IDENSB),ICORE(IDENSBCUB),1,ICORE(IMOBCUB),
     $              ICORE(IPVAL2),ICORE(IIPCNT),ICORE(IAPDX),
     $              ICORE(IAPDY),ICORE(IAPDZ),ICORE(IAPDR),
     $              ICORE(IDENSNMAT),
     $              ICORE(IPCOEFSO),ICORE(IFINMO),ICORE(IFINSALC),
     $              IBIGLOOP,
     $              JBIGLOOP,IDORO,NORBOUTB)
C     
C     Print final results
               CALL CUBEFORM (ICORE(IDENSACUB),ICORE(IDENSBCUB),NPTSZ,
     $              ICORE(IDENSTCUB),ICORE(IMOACUB),ICORE(IMOBCUB),
     $              NBASP,IDORO,IAORBD,IAORBO,ICORE(IMOOUT),NBASPEFF,
     $              NORBOUTA,NORBOUTB,IORBOUT,IDENOUT)
C     
 101        CONTINUE
 100     CONTINUE
 911  CONTINUE
C     
C     If requested, calculate the density difference cube file.
      IF (IDDIFF.EQ.1) THEN
         CALL DIF ()
      END IF
C     
      CALL CRAPSO
C     
      STOP
      END

      
