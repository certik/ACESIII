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
         SUBROUTINE  NLO__GENER_NLO_ORBITALS
     +
     +                    ( IMAX,ZMAX,
     +                      MEMORY,
     +                      NBAS,NATOM,
     +                      MXSHELL,MXNAL,
     +                      BONDSIZE,RYD2HYB,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      ZATOM,XYZCOOR,
     +                      ALPHA,BETA,
     +                      DENSITY,OVERLAP,
     +                      MJUMP,
     +                      SPHERIC,
     +                      SYMMETRY,
     +                      ONLYNAO,ONLYNHO,ONLYNBO,ONLYNLMO,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      ICORE,ZCORE,
     +
     +                              BDSIZE,
     +                              BDATOM,
     +                              COEFFS )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NLO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : Master routine for finding natural localized orbitals
C                within a molecular environment using the first order
C                density matrix and the overlap matrix.
C
C                The density matrix can be of pure spin alpha or beta
C                in case of UHF functions or of mixed spin in case of
C                RHF functions. To determine the separate alpha and
C                beta natural localized orbitals simply run the routine
C                twice with the corresponding alpha and beta first order
C                density matrix.
C
C
C                  Input:
C
C                    IMAX,ZMAX    =  maximum integer + flp memory
C                    MEMORY       =  if true, the code will determine
C                                    just the maximum amount of memory
C                                    needed. Nothing else will be done.
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MXNAL        =  maximum size of atomic l-shell
C                                    space. The atomic l-shell space
C                                    is the total # of contractions for
C                                    an atomic l-shell.
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    should be included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                                    If false, only the Valence NAO
C                                    space will be used.
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A (in increasing order!).
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    ZATOM (A)    =  atomic number for atom A.
C                    XYZCOOR      =  x,y,z-coordinates for all the
C                                    atoms in 3 x NATOM matrix format.
C                    ALPHA        =  if true, the input AO density
C                                    matrix contains the alpha spin
C                                    electron density.
C                    BETA         =  if true, the input AO density
C                                    matrix contains the beta spin
C                                    electron density.
C                    DENSITY      =  full NBAS x NBAS matrix containing
C                                    the expansion coeffs of the density
C                                    in terms of products of the AO
C                                    basis functions. Can contain the
C                                    alpha and/or the beta spin electron
C                                    density (example UHF or RHF), which
C                                    is controlled by the keywords
C                                    ALPHA and BETA. Note that when
C                                    both these keywords are passed as
C                                    false, the routine stops with a
C                                    message.
C                    OVERLAP      =  full NBAS x NBAS overlap matrix.
C                    MJUMP        =  is .true., if the m values in the
C                                    m-space are ordered such that the
C                                    same m values are separated. This
C                                    keyword is necessary because some
C                                    AO basis functions are m-ordered
C                                    differently within each l-shell.
C                    SPHERIC      =  is true, if the l-shells are
C                                    spherical, false if they are
C                                    cartesian.
C                    SYMMETRY     =  if true, the code will attempt
C                                    the production of rotationally
C                                    invariant NBOs and NLMOs.
C                    ONLYx        =  is true, if only x=NAO,NHO,NBO or
C                                    NLMO orbitals are wanted. Note,
C                                    that if any of the ONLYx is true,
C                                    then all those ONLYx preceeding
C                                    that particular ONLYx are forced
C                                    to be false.
C                    MXCHOOSE     =  maximum # of bonds selected to
C                                    be chosen. The maximum is build
C                                    from all # of chosen bonds for all
C                                    bondsizes.
C                    NCHOOSE (B)  =  # of bonds to be chosen for bonds
C                                    of size B. Four cases:
C                                    1) = 9999 => skip search for bonds
C                                       of size B.
C                                    2) = 0 => complete search for all
C                                       possible bonds of size B will
C                                       be performed.
C                                    3) = -n => only n bonds of size B
C                                       will be searched between those
C                                       atomic indices as provided by
C                                       the CHOOSE array.
C                                    4) = +n => same as case 3) but
C                                       followed by a complete search
C                                       for all possible remaining
C                                       bonds of size B.
C                                    Priority level: 1) > 3) > 4) > 2).
C                    CHOOSE       =  element CHOOSE (I,N,B) contains
C                                    the I-th atomic index of the N-th
C                                    chosen bond of size B. The order
C                                    of the atomic indices is arbitrary.
C                    ICORE,ZCORE  =  integer/flp scratch space
C
C
C                  Output:
C
C                    BDSIZE (I)   =  Natural localized orbital bond
C                                    sizes. Indicates the # of atoms
C                                    which form the I-th natural
C                                    localized orbital. 
C                    BDATOM (I,J) =  Atomic index map for the natural
C                                    localized orbitals. J-th atomic
C                                    index forming the I-th natural
C                                    localized orbital.
C                    COEFFS (I,J) =  Natural localized orbitals
C                                    coefficient matrix in AO basis.
C                                    I-th AO coefficient in the J-th
C                                    natural localized orbital.
C                                    
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ALPHA,BETA
         LOGICAL     CYCLERYD,CYCLESYM
         LOGICAL     MEMORY
         LOGICAL     MJUMP
         LOGICAL     ONLYNAO,ONLYNHO,ONLYNBO,ONLYNLMO
         LOGICAL     PLATONIC
         LOGICAL     RYD2HYB
         LOGICAL     SETWPRE
         LOGICAL     SPHERIC
         LOGICAL     SYMMETRY

         INTEGER     BONDSIZE,NBOSIZE,ROTSIZE
         INTEGER     GPOINT
         INTEGER     IMAX,ZMAX
         INTEGER     MXCHOOSE
         INTEGER     MXNCBA,MXNVBA,MXNRBA,MXNHBA,MXNBA,MXCOL
         INTEGER     MXSHELL,MXLSIZE,MXNAL,MX2CEN
         INTEGER     NBAS,NATOM,NBOND
         INTEGER     NCA,NVA,NRA,NHA,NLA,NEA,NYA
         INTEGER     NCYCLE
         INTEGER     NOCC
         INTEGER     NMB,NCB,NVB,NRB,NHB,NLB,NBB,NEB,NAB,NYB
         INTEGER     NTETRA,NCUBE,NOCTA,NDODECA,NICOSA
         INTEGER     UNITID,UNITDB

         INTEGER     IAT2CEN,ICOLMAP,INRYDAL,ILSIZE,INAOTYP,ICSHELL,
     +               IVSHELL,IRSHELL,IHSHELL,IATNCB,IATCIDX,IATCOFF,
     +               IATNVB,IATVIDX,IATVOFF,IATNRB,IATRIDX,IATROFF,
     +               IATNLB,IATLIDX,IATLOFF,IATNEB,IATEIDX,IATEOFF,
     +               IATNYB,IATYIDX,IATYOFF,IATNHB,IATHVAL,IATHIDX,
     +               IATHOFF,IATHCEN,IATORD,ILOCMAP,IVALIDX,IRYDIDX,
     +               INHYB,IBDNCEN,IBDCEN,IBDNBAS,IBDBAS,IBDOCC,INBOBD,
     +               INWSTEP,ISYMNBO,IRING,INRING,IRINGSZ,
     +               IPLATO,IBASBEG,IBASEND,ISYMMAP,ISYMCEN,ISLTYPE,
     +               IANGSYM,IWORK1,IWORK2,IWORK3,IWORK4,ITOP

         INTEGER     ZB,ZH,ZQ,ZW,
     +               ZPA,ZPH,
     +               ZCAL,ZSAL,ZWAL,ZSAH,ZHYB,ZROT,
     +               ZWPRE,ZWNAO,ZTSYM,
     +               ZDOMAT,ZSOMAT,ZBOMAT,ZWBOND,ZWSTAR,ZLOCAL,
     +               ZPHSUB,ZPHDEP,
     +               ZPOPCOR,ZPOPVAL,ZPOPRYD,ZPOPSUM,ZDENORM,
     +               ZANGNHO,ZANGNBO,ZWRYDAT,
     +               ZWORK1,ZWORK2,ZWORK3,ZWORK4,ZWORK5,
     +               ZTOP

         INTEGER     BDSIZE  (1:NBAS)
         INTEGER     ICORE   (1:IMAX)
         INTEGER     NCHOOSE (1:BONDSIZE)
         INTEGER     NSHELLS (1:NATOM)
         INTEGER     ZATOM   (1:NATOM)

         INTEGER     BDATOM  (1:NBAS     ,1:BONDSIZE)
         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM   )
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM   )

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE,1:BONDSIZE)

         DOUBLE PRECISION     DSYMACC,LSYMACC,PSYMACC,QSYMACC
         DOUBLE PRECISION     MAXOCC
         DOUBLE PRECISION     NO2CEN
         DOUBLE PRECISION     ONE
         DOUBLE PRECISION     WBDMIN,WSTMAX,WBDCRT,WSTEP
         DOUBLE PRECISION     WPRERYD,WNAOVAL,WNAORYD,WNBOOCC

         DOUBLE PRECISION     ZCORE   (1:ZMAX)

         DOUBLE PRECISION     COEFFS  (1:NBAS,1:NBAS )
         DOUBLE PRECISION     DENSITY (1:NBAS,1:NBAS )
         DOUBLE PRECISION     OVERLAP (1:NBAS,1:NBAS )
         DOUBLE PRECISION     XYZCOOR (1:3   ,1:NATOM)

         PARAMETER  (ONE     = 1.D0)
         PARAMETER  (ROTSIZE = 3   )
C
C
C------------------------------------------------------------------------
C
C
C             ...open general printout file and diagnostic/debug file.
C
C
         IF (ALPHA.AND.BETA) THEN
             IF (.NOT.MEMORY) THEN
                 CALL  NLO__OPEN_FILES
     +
     +                      ( 'NLO-results-printout',
     +                        'NLO-diagnostic-debug',
     +
     +                                   UNITID,
     +                                   UNITDB )
     +
     +
             END IF
         ELSE IF (ALPHA) THEN
             IF (.NOT.MEMORY) THEN
                 CALL  NLO__OPEN_FILES
     +
     +                      ( 'NLO-results-printout-alpha',
     +                        'NLO-diagnostic-debug-alpha',
     +
     +                                   UNITID,
     +                                   UNITDB )
     +
     +
             END IF
         ELSE IF (BETA) THEN
             IF (.NOT.MEMORY) THEN
                 CALL  NLO__OPEN_FILES
     +
     +                      ( 'NLO-results-printout-beta',
     +                        'NLO-diagnostic-debug-beta',
     +
     +                                   UNITID,
     +                                   UNITDB )
     +
     +
             END IF
         ELSE
             WRITE (*,*) ' ALPHA and BETA set to false! '
             WRITE (*,*) ' Unable to proceed! '
             WRITE (*,*) ' nlo__gener_nlo_orbitals '
             STOP
         END IF
C
C
C             ...preliminary steps.
C
C
         IF (ONLYNLMO) THEN
             ONLYNAO = .FALSE.
             ONLYNHO = .FALSE.
             ONLYNBO = .FALSE.
         END IF

         IF (ONLYNBO) THEN
             ONLYNAO = .FALSE.
             ONLYNHO = .FALSE.
         END IF

         IF (ONLYNHO) THEN
             ONLYNAO = .FALSE.
         END IF

         IF (ALPHA.AND.BETA) THEN
             MAXOCC = ONE + ONE
         ELSE
             MAXOCC = ONE
         END IF

         MX2CEN  = NATOM * (NATOM + 1) / 2

         IBASBEG = 1
         IBASEND = IBASBEG + NATOM
         IATORD  = IBASEND + NATOM
         ILSIZE  = IATORD  + NATOM
         INRING  = ILSIZE  + (MXSHELL + 1)
         IRINGSZ = INRING  + NATOM
         IRING   = IRINGSZ + NATOM
         IPLATO  = IRING   + NATOM * NATOM
         ISYMMAP = IPLATO  + NATOM * 5
         IAT2CEN = ISYMMAP + NATOM * NATOM
         IWORK1  = IAT2CEN + 2*MX2CEN
         ITOP    = IWORK1  + MAX (NBAS,NATOM,MX2CEN)

         ZDENORM = 1
         ZDOMAT  = ZDENORM + NBAS
         ZSOMAT  = ZDOMAT  + NATOM * NATOM
         ZBOMAT  = ZSOMAT  + NATOM * NATOM
         ZLOCAL  = ZBOMAT  + NATOM * NATOM
         ZWORK1  = ZLOCAL  + NATOM * NBAS
         ZWORK2  = ZWORK1  + NBAS
         ZTOP    = ZWORK2  + NBAS * NBAS

         IF (MEMORY) THEN
             IMAX = ITOP
             ZMAX = ZTOP
         ELSE
             CALL  NLO__NORMALIZE_OVERLAP_DENSITY
     +
     +                  ( NBAS,
     +                    ZCORE (ZWORK1),
     +                    ICORE (IWORK1),
     +
     +                              ZCORE (ZDENORM),
     +                              DENSITY,OVERLAP )
     +
     +
             CALL  NLO__INITIALIZE_RUN
     +
     +                  ( NBAS,NATOM,
     +                    MXSHELL,
     +                    MAXOCC,BONDSIZE,
     +                    NSHELLS,SHELLS,NBASAL,
     +                    SPHERIC,
     +                    MXCHOOSE,NCHOOSE,CHOOSE,
     +                    ICORE (IWORK1),
     +                    XYZCOOR,ZATOM,
     +                    OVERLAP,
     +                    ZCORE (ZWORK1),
     +                    ZCORE (ZWORK2),
     +
     +                              NO2CEN,
     +                              WPRERYD,WNAOVAL,WNAORYD,WNBOOCC,
     +                              WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                              DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                              MXNBA,MXLSIZE,
     +                              ICORE (ILSIZE),
     +                              ICORE (IBASBEG),
     +                              ICORE (IBASEND),
     +                              ICORE (ISYMMAP),
     +                              BDATOM,
     +                              GPOINT,
     +                              ZCORE (ZDOMAT),
     +                              ZCORE (ZSOMAT),
     +                              ZCORE (ZBOMAT),
     +                              PLATONIC,
     +                              DENSITY )
     +
     +
             CALL  NLO__PRINT_MOLECULE_INFO
     +
     +                  ( UNITID,
     +                    NATOM,
     +                    XYZCOOR,ZATOM,
     +                    GPOINT,
     +                    ZCORE (ZDOMAT),
     +                    ZCORE (ZSOMAT),
     +                    ZCORE (ZBOMAT) )
     +
     +
             CALL  NLO__ANALYZE_BOND_ORDER_MATRIX
     +
     +                  ( NATOM,
     +                    MX2CEN,
     +                    ICORE (IAT2CEN),
     +                    NO2CEN,
     +                    ICORE (IWORK1),
     +                    ZCORE (ZWORK1),
     +
     +                              ICORE (IATORD),
     +                              ZCORE (ZBOMAT) )
     +
     +
         END IF
C
C
C             ...the NAO generation section follows.
C
C
         IANGSYM = ISYMMAP + NATOM * NATOM
         ICOLMAP = IANGSYM + NBAS
         INRYDAL = ICOLMAP + NBAS
         INAOTYP = INRYDAL + NATOM * (MXSHELL + 1)
         ICSHELL = INAOTYP + NBAS
         IVSHELL = ICSHELL + NBAS
         IRSHELL = IVSHELL + NBAS
         IATNCB  = IRSHELL + NBAS
         IATCIDX = IATNCB  + NATOM
         IATCOFF = IATCIDX + NATOM
         IATNVB  = IATCOFF + NATOM
         IATVIDX = IATNVB  + NATOM
         IATVOFF = IATVIDX + NATOM
         IATNRB  = IATVOFF + NATOM
         IATRIDX = IATNRB  + NATOM
         IATROFF = IATRIDX + NATOM
         IWORK1  = IATROFF + NATOM
         IWORK2  = IWORK1  + NATOM
         IWORK3  = IWORK2  + NATOM
         ITOP    = IWORK3  + NBAS + NBAS

         ZWRYDAT = ZLOCAL  + NATOM * NBAS
         ZWNAO   = ZWRYDAT + NATOM
         ZWPRE   = ZWNAO   + NBAS
         ZSAL    = ZWPRE   + NBAS
         ZCAL    = ZSAL    + MXNAL * MXNAL
         ZWAL    = ZCAL    + MXNAL * MXNAL
         ZTSYM   = ZWAL    + MXNAL
         ZPA     = ZTSYM   + MXLSIZE * MXLSIZE
         ZPOPCOR = ZPA     + MXNBA * MXNBA
         ZPOPVAL = ZPOPCOR + NATOM
         ZPOPRYD = ZPOPVAL + NATOM
         ZPOPSUM = ZPOPRYD + NATOM
         ZWORK1  = ZPOPSUM + NATOM
         ZWORK2  = ZWORK1  + NBAS + NBAS
         ZTOP    = ZWORK2  + NBAS * NBAS

         IF (MEMORY) THEN
             IMAX = MAX (IMAX,ITOP)
             ZMAX = MAX (ZMAX,ZTOP)
         ELSE
             NCYCLE = 0
             SETWPRE = .TRUE.

 1000        NCYCLE = NCYCLE + 1
             WRITE (*,*) ' pre-NAO/NAO cycle # ',NCYCLE

             CALL  NLO__GENER_PRE_NAO_ORBITALS
     +
     +                  ( NBAS,NATOM,
     +                    MXSHELL,MXNAL,
     +                    NSHELLS,SHELLS,NBASAL,
     +                    ICORE (ILSIZE),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    ICORE (ISYMMAP),
     +                    ICORE (IWORK1),
     +                    ICORE (IWORK2),
     +                    DENSITY,OVERLAP,
     +                    ZCORE (ZSAL),
     +                    ZCORE (ZCAL),
     +                    ZCORE (ZWAL),
     +                    SETWPRE,
     +                    WPRERYD,
     +                    ZCORE (ZWRYDAT),
     +                    MJUMP,
     +
     +                             NMB,NRB,
     +                             ICORE (ICOLMAP),
     +                             ICORE (INRYDAL),
     +                             ZCORE (ZWPRE),
     +                             COEFFS )
     +
     +
             CALL  NLO__GENER_NAO_ORBITALS
     +
     +                  ( NBAS,NATOM,
     +                    MXSHELL,MXLSIZE,MXNAL,MXNBA,
     +                    NSHELLS,SHELLS,NBASAL,
     +                    DENSITY,
     +                    ZCORE (ZPA),
     +                    OVERLAP,
     +                    ZCORE (ZSAL),
     +                    ZCORE (ZCAL),
     +                    ZCORE (ZWAL),
     +                    ZCORE (ZTSYM),
     +                    MJUMP,
     +                    NMB,NRB,
     +                    ICORE (ILSIZE),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    PSYMACC,
     +                    ICORE (ISYMMAP),
     +                    ICORE (IWORK1),
     +                    ICORE (IWORK2),
     +                    ICORE (ICOLMAP),
     +                    ICORE (INRYDAL),
     +                    ICORE (IWORK3),
     +                    ZCORE (ZWORK1),
     +                    ZCORE (ZWORK2),
     +                    ZCORE (ZWPRE),
     +
     +                             BDSIZE,
     +                             BDATOM (1,1),
     +                             ICORE (IANGSYM),
     +                             ZCORE (ZWNAO),
     +                             COEFFS )
     +
     +
             CALL  NLO__ANALYZE_NAO_ORBITALS
     +
     +                  ( NBAS,NATOM,
     +                    MXSHELL,MXNAL,
     +                    MAXOCC,
     +                    NSHELLS,SHELLS,NBASAL,
     +                    ZATOM,
     +                    NMB,
     +                    ICORE (ILSIZE),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    PSYMACC,
     +                    WNAOVAL,WNAORYD,
     +                    ZCORE (ZWPRE),
     +                    ZCORE (ZWNAO),
     +                    ZCORE (ZWORK1),
     +
     +                             NCB,NVB,NRB,
     +                             ZCORE (ZWRYDAT),
     +                             ZCORE (ZPOPCOR),
     +                             ZCORE (ZPOPVAL),
     +                             ZCORE (ZPOPRYD),
     +                             ZCORE (ZPOPSUM),
     +                             ICORE (INAOTYP),
     +                             ICORE (ICOLMAP),
     +                             ICORE (ISYMMAP),
     +                             ICORE (ICSHELL),
     +                             ICORE (IVSHELL),
     +                             ICORE (IRSHELL),
     +                             NCA,
     +                             MXNCBA,
     +                             ICORE (IATNCB),
     +                             ICORE (IATCIDX),
     +                             ICORE (IATCOFF),
     +                             NVA,
     +                             MXNVBA,
     +                             ICORE (IATNVB),
     +                             ICORE (IATVIDX),
     +                             ICORE (IATVOFF),
     +                             NRA,
     +                             MXNRBA,
     +                             ICORE (IATNRB),
     +                             ICORE (IATRIDX),
     +                             ICORE (IATROFF),
     +                             CYCLERYD,
     +                             CYCLESYM )
     +
     +
             SETWPRE = CYCLESYM
             IF (CYCLERYD .OR. CYCLESYM) GOTO 1000

         END IF
C
C
C             ...symmetry analysis based on the final NAOs obtained.
C
C
         IWORK1  = IATROFF + NATOM
         IWORK2  = IWORK1  + NATOM
         IWORK3  = IWORK2  + NATOM
         IWORK4  = IWORK3  + NATOM
         ITOP    = IWORK4  + NATOM

         ZWORK1  = ZPOPSUM + NATOM
         ZWORK2  = ZWORK1  + NATOM
         ZWORK3  = ZWORK2  + NATOM
         ZWORK4  = ZWORK3  + NATOM
         ZWORK5  = ZWORK4  + NATOM
         ZTOP    = ZWORK5  + NATOM

         IF (MEMORY) THEN
             IMAX = MAX (IMAX,ITOP)
             ZMAX = MAX (ZMAX,ZTOP)
         ELSE
             CALL  NLO__SYMMETRY_RELATED_CENTERS
     +
     +                  ( NATOM,
     +                    ICORE (IWORK1),
     +                    ICORE (IWORK2),
     +                    ICORE (IWORK3),
     +                    ICORE (ISYMMAP),
     +                    ZCORE (ZDOMAT),
     +                    DSYMACC,
     +                    XYZCOOR,
     +                    ICORE (IWORK4),
     +                    ZCORE (ZWORK1),
     +                    ZCORE (ZWORK2),
     +                    ZCORE (ZWORK3),
     +                    ZCORE (ZWORK4),
     +                    ZCORE (ZWORK5),
     +
     +                             PLATONIC,
     +                             ICORE (IRING),
     +                             ICORE (INRING),
     +                             ICORE (IRINGSZ),
     +                             NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                             ICORE (IPLATO) )
     +
     +
             CALL  NLO__PRINT_SYMMETRY_INFO
     +
     +                  ( UNITID,
     +                    NATOM,
     +                    ZATOM,
     +                    ICORE (ISYMMAP),
     +                    ICORE (IRING),
     +                    ICORE (INRING),
     +                    ICORE (IRINGSZ),
     +                    ICORE (IPLATO) )
     +
     +
         END IF
C
C
C             ...generate the square root of the overlap matrix
C                needed for NAO,NBO and NLMO locality analysis.
C                The overlap matrix is not needed any more and
C                is thus overwritten by its square root.
C
C
         ZWORK1  = ZPOPSUM + NATOM
         ZWORK2  = ZWORK1  + NBAS
         ZTOP    = ZWORK2  + NBAS * NBAS

         IF (MEMORY) THEN
             ZMAX = MAX (ZMAX,ZTOP)
             IF (ONLYNAO) THEN
                 RETURN
             END IF
         ELSE
             CALL  NLO__GENER_SQROOT_OVERLAP
     +
     +                  ( NBAS,
     +                    ZCORE (ZWORK1),
     +                    ZCORE (ZWORK2),
     +
     +                             OVERLAP )
     +
     +
             CALL  NLO__DETERMINE_ORBITAL_LOCALITY
     +
     +                  ( NBAS,NATOM,
     +                    ICORE (IWORK1),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    COEFFS,
     +                    OVERLAP,
     +
     +                             ZCORE (ZLOCAL) )
     +
     +
             CALL  NLO__PRINT_NAO_RESULTS
     +
     +                  ( UNITID,
     +                    NBAS,NATOM,
     +                    MXSHELL,
     +                    NSHELLS,SHELLS,NBASAL,
     +                    ZATOM,
     +                    ICORE (ILSIZE),
     +                    ZCORE (ZPOPCOR),
     +                    ZCORE (ZPOPVAL),
     +                    ZCORE (ZPOPRYD),
     +                    ZCORE (ZPOPSUM),
     +                    ICORE (INAOTYP),
     +                    ICORE (IANGSYM),
     +                    ZCORE (ZLOCAL),
     +                    ZCORE (ZWNAO) )
     +
     +
             IF (ONLYNAO) THEN
                 CALL  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +                      ( NBAS,
     +                        ZCORE (ZDENORM),
     +
     +                                 COEFFS )
     +
     +
                 RETURN
             END IF

         END IF
C
C
C             ...enter the NHO generation.
C
C
         IHSHELL = IATROFF + NATOM
         IATNHB  = IHSHELL + NBAS
         IATHVAL = IATNHB  + NATOM
         IATHIDX = IATHVAL + NATOM
         IATHOFF = IATHIDX + NATOM
         ILOCMAP = IATHOFF + NATOM
         IVALIDX = ILOCMAP + NBAS
         IRYDIDX = IVALIDX + NATOM
         IWORK1  = IRYDIDX + NATOM
         ITOP    = IWORK1  + NBAS

         IF (MEMORY) THEN
             IMAX = MAX (IMAX,ITOP)
         ELSE
             CALL  NLO__CHARACTERIZE_NHO_SPACE
     +
     +                  ( NBAS,NATOM,
     +                    NVB,NVA,
     +                    MXNVBA,
     +                    ICORE (IVSHELL),
     +                    ICORE (IATNVB),
     +                    ICORE (IATVIDX),
     +                    ICORE (IATVOFF),
     +                    RYD2HYB,
     +                    ICORE (ILOCMAP),
     +                    ICORE (IVALIDX),
     +                    ICORE (IRYDIDX),
     +                    ICORE (IWORK1),
     +
     +                             ICORE (IATORD),
     +                             NHB,NCB,NRB,
     +                             ICORE (ICOLMAP),
     +                             ICORE (IHSHELL),
     +                             ICORE (ICSHELL),
     +                             ICORE (IRSHELL),
     +                             NHA,
     +                             MXNHBA,
     +                             ICORE (IATNHB),
     +                             ICORE (IATHVAL),
     +                             ICORE (IATHIDX),
     +                             ICORE (IATHOFF),
     +                             NCA,
     +                             MXNCBA,
     +                             ICORE (IATNCB),
     +                             ICORE (IATCIDX),
     +                             ICORE (IATCOFF),
     +                             NRA,
     +                             MXNRBA,
     +                             ICORE (IATNRB),
     +                             ICORE (IATRIDX),
     +                             ICORE (IATROFF),
     +                             MXNBA,MXCOL )
     +
     +
         END IF

         IF (SYMMETRY) THEN
             NBOSIZE = MAX (BONDSIZE,ROTSIZE)
         ELSE
             NBOSIZE = BONDSIZE
         END IF

         IATHCEN = IATHOFF + NATOM
         INHYB   = IATHCEN + NATOM
         ISLTYPE = INHYB   + NATOM
         IBDNBAS = ISLTYPE + NBAS
         IBDNCEN = IBDNBAS + NBAS
         IBDOCC  = IBDNCEN + NBAS
         IBDCEN  = IBDOCC  + NBAS
         IBDBAS  = IBDCEN  + NBOSIZE * NBAS
         INWSTEP = IBDBAS  + NBOSIZE * NBAS
         IWORK1  = INWSTEP + BONDSIZE
         ITOP    = IWORK1  + 2 * (NBAS + NHB + NHA)

         ZW      = ZWNAO
         ZH      = ZW + NBAS
         ZANGNHO = ZH + MXNBA * NBAS
         ZWBOND  = ZANGNHO + NBAS * (MXSHELL + 1)
         ZWSTAR  = ZWBOND + BONDSIZE
         ZSAH    = ZWSTAR + BONDSIZE
         ZPH     = ZSAH + MXNHBA * MXNHBA
         ZPHSUB  = ZPH + NHB * NHB
         ZPHDEP  = ZPHSUB + NHB * NHB
         ZWORK1  = ZPHDEP + NHB * NHB
         ZWORK2  = ZWORK1 + 3 * NBAS
         ZTOP    = ZWORK2 + NBAS * MXCOL

         IF (MEMORY) THEN
             IMAX = MAX (IMAX,ITOP)
             ZMAX = MAX (ZMAX,ZTOP)
         ELSE
             CALL  NLO__GENER_NHO_ORBITALS
     +
     +                  ( NBAS,NATOM,
     +                    MXNHBA,MXNBA,MXSHELL,MXCOL,
     +                    BONDSIZE,NBOSIZE,
     +                    MAXOCC,
     +                    NHB,NCB,NRB,
     +                    NHA,NCA,NRA,
     +                    ICORE (IATNHB),
     +                    ICORE (IATHVAL),
     +                    ICORE (IATHOFF),
     +                    ICORE (IATNCB),
     +                    ICORE (IATCIDX),
     +                    ICORE (IATCOFF),
     +                    ICORE (IATNRB),
     +                    ICORE (IATRIDX),
     +                    ICORE (IATROFF),
     +                    ICORE (IATORD),
     +                    ICORE (IATHCEN),
     +                    ICORE (IHSHELL),
     +                    ICORE (ICSHELL),
     +                    ICORE (IRSHELL),
     +                    RYD2HYB,
     +                    ICORE (ICOLMAP),
     +                    ICORE (INHYB),
     +                    ICORE (INWSTEP),
     +                    DENSITY,
     +                    ZCORE (ZSAH),
     +                    ZCORE (ZPH),
     +                    ZCORE (ZPHSUB),
     +                    ZCORE (ZPHDEP),
     +                    ZCORE (ZBOMAT),
     +                    WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                    MXCHOOSE,NCHOOSE,CHOOSE,
     +                    ICORE (IWORK1),
     +                    ZCORE (ZWORK1),
     +                    ZCORE (ZWORK2),
     +
     +                             NBOND,
     +                             ZCORE (ZWBOND),
     +                             ZCORE (ZWSTAR),
     +                             ICORE (IBDNCEN),
     +                             ICORE (IBDCEN),
     +                             ICORE (IBDNBAS),
     +                             ICORE (IBDBAS),
     +                             ICORE (IBDOCC),
     +                             BDATOM (1,1),
     +                             ZCORE (ZANGNHO),
     +                             ZCORE (ZH),
     +                             ZCORE (ZW),
     +                             COEFFS )
     +
     +
             CALL  NLO__PRINT_NHO_RESULTS
     +
     +                  ( UNITID,
     +                    NBAS,NATOM,NBOND,
     +                    NBOSIZE,
     +                    MXNBA,MXSHELL,
     +                    NHB,NCB,NRB,
     +                    NHA,NCA,NRA,
     +                    ZATOM,
     +                    RYD2HYB,
     +                    ICORE (IATNHB),
     +                    ICORE (IATHVAL),
     +                    ICORE (IATHIDX),
     +                    ICORE (IATHOFF),
     +                    ICORE (IATNCB),
     +                    ICORE (IATCIDX),
     +                    ICORE (IATCOFF),
     +                    ICORE (IATNRB),
     +                    ICORE (IATRIDX),
     +                    ICORE (IATROFF),
     +                    ICORE (IHSHELL),
     +                    ICORE (ICSHELL),
     +                    ICORE (IRSHELL),
     +                    ICORE (IBDNCEN),
     +                    ICORE (IBDCEN),
     +                    ICORE (IBDNBAS),
     +                    ICORE (IBDBAS),
     +                    ICORE (IBDOCC),
     +                    ICORE (ICOLMAP),
     +                    ICORE (IWORK1),
     +                    ZCORE (ZANGNHO),
     +                    ZCORE (ZH) )
     +
     +
         END IF

         IF (ONLYNHO) THEN
             RETURN
         END IF
C
C
C             ...enter the NBO generation.
C
C
         INBOBD  = IBDBAS  + NBOSIZE * NBAS
         ISYMNBO = INBOBD  + NBAS
         ISYMCEN = ISYMNBO + NBAS
         IATNLB  = ISYMCEN + NATOM
         IATLIDX = IATNLB  + NATOM
         IATLOFF = IATLIDX + NATOM
         IATNEB  = IATLOFF + NATOM
         IATEIDX = IATNEB  + NATOM
         IATEOFF = IATEIDX + NATOM
         IATNYB  = IATEOFF + NATOM
         IATYIDX = IATNYB  + NATOM
         IATYOFF = IATYIDX + NATOM
         IWORK1  = IATYOFF + NATOM
         ITOP    = IWORK1  + MAX (3*NBAS,NATOM)

         ZQ      = ZWSTAR + BONDSIZE
         ZB      = ZQ + NBAS
         ZANGNBO = ZB + NBOSIZE * NBAS
         ZPH     = ZANGNBO + NBAS * (MXSHELL + 1)
         ZHYB    = ZPH + NHB * NHB
         ZROT    = ZHYB + ROTSIZE * ROTSIZE
         ZWORK1  = ZROT + ROTSIZE * ROTSIZE
         ZWORK2  = ZWORK1 + NBAS + NBAS
         ZTOP    = ZWORK2 + NBAS * NBAS

         IF (MEMORY) THEN
             IMAX = MAX (IMAX,ITOP)
             ZMAX = MAX (ZMAX,ZTOP)
             IF (ONLYNBO) THEN
                 RETURN
             END IF
         ELSE
             CALL  NLO__GENER_NBO_ORBITALS
     +
     +                  ( NBAS,NATOM,NBOND,
     +                    BONDSIZE,NBOSIZE,
     +                    MXNBA,MXSHELL,
     +                    NHB,NCB,NRB,
     +                    NHA,NCA,NRA,
     +                    ICORE (IATNHB),
     +                    ICORE (IATHOFF),
     +                    ICORE (IATNCB),
     +                    ICORE (IATCIDX),
     +                    ICORE (IATCOFF),
     +                    ICORE (IATNRB),
     +                    ICORE (IATRIDX),
     +                    ICORE (IATROFF),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    ICORE (ICOLMAP),
     +                    ZCORE (ZWBOND),
     +                    ZCORE (ZWSTAR),
     +                    WNAORYD,WNBOOCC,
     +                    ICORE (IBDNCEN),
     +                    ICORE (IBDCEN),
     +                    ICORE (IBDNBAS),
     +                    ICORE (IBDBAS),
     +                    ICORE (IBDOCC),
     +                    ZCORE (ZANGNHO),
     +                    DENSITY,
     +                    ZCORE (ZH),
     +                    ZCORE (ZPH),
     +                    ICORE (IWORK1),
     +                    ZCORE (ZWORK1),
     +                    ZCORE (ZWORK2),
     +
     +                             NLA,NEA,NYA,
     +                             NLB,NBB,NEB,NAB,NYB,
     +                             ICORE (IATNLB),
     +                             ICORE (IATLIDX),
     +                             ICORE (IATLOFF),
     +                             ICORE (IATNEB),
     +                             ICORE (IATEIDX),
     +                             ICORE (IATEOFF),
     +                             ICORE (IATNYB),
     +                             ICORE (IATYIDX),
     +                             ICORE (IATYOFF),
     +                             ICORE (INBOBD),
     +                             ZCORE (ZANGNBO),
     +                             BDSIZE,BDATOM,
     +                             ZCORE (ZB),
     +                             ZCORE (ZW),
     +                             ZCORE (ZQ),
     +                             COEFFS )
     +
     +
             IF (SYMMETRY) THEN

                 CALL  NLO__SYMMETRIZE_NBO_ORBITALS
     +
     +                      ( NBAS,NATOM,NBOND,
     +                        NBOSIZE,ROTSIZE,
     +                        MXNBA,MXSHELL,
     +                        NCB,NLB,NBB,NEB,NAB,NYB,NRB,
     +                        NCA,NLA,NEA,NYA,NRA,
     +                        ICORE (IATNCB),
     +                        ICORE (IATCIDX),
     +                        ICORE (IATCOFF),
     +                        ICORE (IATNLB),
     +                        ICORE (IATLIDX),
     +                        ICORE (IATLOFF),
     +                        ICORE (IATNEB),
     +                        ICORE (IATEIDX),
     +                        ICORE (IATEOFF),
     +                        ICORE (IATNYB),
     +                        ICORE (IATYIDX),
     +                        ICORE (IATYOFF),
     +                        ICORE (IATNRB),
     +                        ICORE (IATRIDX),
     +                        ICORE (IATROFF),
     +                        ICORE (IBASBEG),
     +                        ICORE (IBASEND),
     +                        DSYMACC,LSYMACC,
     +                        PSYMACC,QSYMACC,
     +                        ICORE (ISYMMAP),
     +                        ICORE (ISYMNBO),
     +                        ICORE (ISYMCEN),
     +                        PLATONIC,
     +                        ICORE (IRING),
     +                        ICORE (INRING),
     +                        ICORE (IRINGSZ),
     +                        NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                        ICORE (IPLATO),
     +                        ICORE (IBDNCEN),
     +                        ICORE (IBDCEN),
     +                        ZCORE (ZANGNHO),
     +                        ZCORE (ZHYB),
     +                        ZCORE (ZROT),
     +                        DENSITY,
     +                        ZCORE (ZDOMAT),
     +                        OVERLAP,
     +                        ZCORE (ZLOCAL),
     +                        ICORE (ISLTYPE),
     +                        ICORE (IWORK1),
     +                        ZCORE (ZWORK1),
     +                        ZCORE (ZWORK2),
     +
     +                                 ICORE (INBOBD),
     +                                 ICORE (IBDNBAS),
     +                                 ICORE (IBDBAS),
     +                                 ZCORE (ZANGNBO),
     +                                 ZCORE (ZB),
     +                                 ZCORE (ZW),
     +                                 ZCORE (ZQ),
     +                                 COEFFS )
     +
     +
             END IF

             CALL  NLO__DETERMINE_ORBITAL_LOCALITY
     +
     +                  ( NBAS,NATOM,
     +                    ICORE (IWORK1),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    COEFFS,
     +                    OVERLAP,
     +
     +                             ZCORE (ZLOCAL) )
     +
     +
             CALL  NLO__PRINT_NBO_RESULTS
     +
     +                  ( UNITID,
     +                    NBAS,NATOM,
     +                    NBOSIZE,
     +                    MXSHELL,
     +                    NCA,NLA,NEA,NYA,NRA,
     +                    NCB,NLB,NBB,NEB,NAB,NYB,NRB,
     +                    ZATOM,
     +                    ICORE (IATNCB),
     +                    ICORE (IATCIDX),
     +                    ICORE (IATCOFF),
     +                    ICORE (IATNLB),
     +                    ICORE (IATLIDX),
     +                    ICORE (IATLOFF),
     +                    ICORE (IATNEB),
     +                    ICORE (IATEIDX),
     +                    ICORE (IATEOFF),
     +                    ICORE (IATNYB),
     +                    ICORE (IATYIDX),
     +                    ICORE (IATYOFF),
     +                    ICORE (IATNRB),
     +                    ICORE (IATRIDX),
     +                    ICORE (IATROFF),
     +                    ICORE (IBDNCEN),
     +                    ICORE (IBDCEN),
     +                    ICORE (IBDNBAS),
     +                    ICORE (IBDBAS),
     +                    ICORE (INBOBD),
     +                    ZCORE (ZLOCAL),
     +                    ZCORE (ZANGNBO),
     +                    ZCORE (ZW),
     +                    ZCORE (ZQ),
     +                    ZCORE (ZB) )
     +
     +
             IF (ONLYNBO) THEN
                 CALL  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +                      ( NBAS,
     +                        ZCORE (ZDENORM),
     +
     +                                 COEFFS )
     +
     +
                 RETURN
             END IF
         END IF
C
C
C             ...enter the NLMO generation.
C
C
         IWORK1  = IBASEND + NATOM
         ITOP    = IWORK1  + MAX (NBAS,NATOM)

         ZW      = ZDENORM + NBAS
         ZLOCAL  = ZW + NBAS
         ZWORK1  = ZLOCAL + NATOM * NBAS
         ZTOP    = ZWORK1 + NBAS * NBAS

         IF (MEMORY) THEN
             IMAX = MAX (IMAX,ITOP)
             ZMAX = MAX (ZMAX,ZTOP)
         ELSE
             NOCC = NCB + NLB + NBB

             CALL  NLO__GENER_NLMO_ORBITALS
     +
     +                  ( NBAS,NOCC,
     +                    MAXOCC,
     +                    ICORE (IWORK1),
     +                    PSYMACC,QSYMACC,
     +                    DENSITY,
     +                    ZCORE (ZWORK1),
     +
     +                             ZCORE (ZW),
     +                             COEFFS )
     +
     +
             CALL  NLO__DETERMINE_ORBITAL_LOCALITY
     +
     +                  ( NBAS,NATOM,
     +                    ICORE (IWORK1),
     +                    ICORE (IBASBEG),
     +                    ICORE (IBASEND),
     +                    COEFFS,
     +                    OVERLAP,
     +
     +                             ZCORE (ZLOCAL) )
     +
     +
             CALL  NLO__PRINT_NLMO_RESULTS
     +
     +                  ( UNITID,
     +                    NBAS,NATOM,
     +                    ZCORE (ZLOCAL) )
     +
     +
             CALL  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +                  ( NBAS,
     +                    ZCORE (ZDENORM),
     +
     +                             COEFFS )
     +
     +
         END IF
C
C
C             ...close general printout file and diagnostic/debug file.
C
C
         IF (.NOT.MEMORY) THEN
             CALL  NLO__CLOSE_FILES
     +
     +                  ( UNITID,UNITDB )
     +
     +
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
