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
         SUBROUTINE  NLO__GENER_NHO_ORBITALS
     +
     +                    ( NBAS,NATOM,
     +                      MXNHBA,MXNBA,MXSHELL,MXCOL,
     +                      BONDSIZE,NBOSIZE,
     +                      MAXOCC,
     +                      NHB,NCB,NRB,
     +                      NHA,NCA,NRA,
     +                      ATNHB,ATHVAL,ATHOFF,
     +                      ATNCB,ATCIDX,ATCOFF,
     +                      ATNRB,ATRIDX,ATROFF,
     +                      ATORD,ATHCEN,
     +                      HSHELL,CSHELL,RSHELL,
     +                      RYD2HYB,
     +                      COLMAP,
     +                      NHYB,
     +                      NWSTEP,
     +                      P,
     +                      SAH,
     +                      PH,PHSUB,PHDEP,
     +                      BOMAT,
     +                      WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      IVEC,XVEC,
     +                      XMAT,
     +
     +                             NBOND,
     +                             WBOND,WSTAR,
     +                             BDNCEN,
     +                             BDCEN,
     +                             BDNBAS,
     +                             BDBAS,
     +                             BDOCC,
     +                             BDATOM,
     +                             ANGNHO,
     +                             H,W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NHO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the sets of natural hybrid
C                orbitals (NHO's) as a set of coefficients in terms
C                of the NHB Hybrid NAOs. The procedure by which all
C                NHB Hybrid NAOs are transformed into NHOs is as
C                follows:
C
C
C                   1) Calculate the NHB x NHB occupation matrix PH
C                      in the Hybrid NAO basis. This matrix will
C                      have the following structure:
C
C
C                             A     B     C     D
C                           -----------------------
C                          |     |     |     |     |
C                        A |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C                          |     |     |     |     |
C                        B |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C                          |     |     |     |     |
C                        C |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C                          |     |     |     |     |
C                        D |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C
C
C                      where A,B,C,D denote atoms.
C
C                   2) Perform the following loop, until the sets of
C                      NHOs found exhaust the NHB space:
C
C                          # of centers NHCEN = 1
C
C                          loop over all NHCEN atomic combinations
C
C                             for each NHCEN atomic combination:
C
C                               i) set up submatrix PHSUB of PH
C                                  by taking appropriate subblocks
C                                  of PH corresponding to the NHCEN
C                                  atomic combination
C
C                              ii) Diagonalize PHSUB and extract
C                                  those eigenfunctions which have
C                                  eigenvalues > Threshold value
C
C                             iii) Deplete PHSUB from those
C                                  eigenfunctions
C
C                              iv) Replace PHSUB section of PH
C                                  by PHSUB matrix
C
C                             if # of NHOs complete, exit
C
C                          continue
C
C                          remove all atomic sites which have by
C                          now a complete set of NHOs.
C
C                          if # of NHOs incomplete, NHCEN = NHCEN + 1
C                          and return to loop beginning
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    MXSHELL      =  largest l-shell value
C                    MXCOL        =  maximum # of NBAS-sized columns
C                                    that are needed for the flp
C                                    scratch array XMAT (see below).
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    NxB          =  total # of Hybrid, Core and
C                                    Rydberg NAOs (x=H,C,R)
C                    NxA          =  total # of hybrid, core and
C                                    Rydberg atoms (x=H,C,R)
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid NAOs
C                                    on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    ATNCB (A)    =  # of Core NAOs on core atom A.
C                    ATCIDX (A)   =  atomic index for core atom A.
C                    ATCOFF (A)   =  index offset for Core NAOs for
C                                    core atom A. This index is equal
C                                    to the total number of Core NAOs
C                                    on all core atoms preceeding
C                                    core atom A.
C                    ATNRB (A)    =  # of Rydberg NAOs on Rydberg
C                                    atom A.
C                    ATRIDX (A)   =  atomic index for Rydberg atom A.
C                    ATROFF (A)   =  index offset for Rydberg NAOs for
C                                    Rydberg atom A. This index is
C                                    equal to the total number of
C                                    Rydberg NAOs on all Rydberg atoms
C                                    preceeding Rydberg atom A.
C                    ATORD (I)    =  I-th atomic index of optimum
C                                    atomic index ordering array of
C                                    all hybrid atoms to be used for
C                                    bond/antibond formation search.
C                    ATHCEN (I)   =  will contain hybrid atomic labels
C                                    (indices) for a subset of all NHA
C                                    atomic hybrid centers
C                    HSHELL (I)   =  l-shell type for I-th Hybrid NAO.
C                    CSHELL (I)   =  l-shell type for I-th Core NAO.
C                    RSHELL (I)   =  l-shell type for I-th Rydberg NAO.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    should be included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                                    If false, only the Valence NAO
C                                    space will be used.
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of I-th NAO, based
C                                    on atomic order, in the NHB/NCB/NRB
C                                    order.
C                    NHYB (A)     =  will contain total # of pre-NHOs
C                                    on each hybrid atom A at each stage
C                                    of the calculation.
C                    NWSTEP (x)   =  will contain the # of steps for
C                                    decreasing/increasing the weight
C                                    acceptance criterion for
C                                    bond/antibond formation for bonds
C                                    over x centers.
C                    P            =  full NBAS x NBAS occupation matrix
C                                    in original AO basis.
C                    SAH          =  will contain the atomic hybrid
C                                    overlap matrices between pre-NHOs.
C                    PH           =  will contain the (possibly depleted)
C                                    NHB x NHB hybrid occupation matrix
C                                    in NAO basis.
C                    PHSUB        =  will contain the submatrices of
C                                    the occupation matrix PH in NAO
C                                    basis for finding the NHOs.
C                    PHDEP        =  will be used for accumulation of
C                                    the depleted occupation matrix PH.
C                    BOMAT (A,B)  =  simplified bond order matrix
C                                    containing info if atoms A and B
C                                    are considered to be bonded or not
C                    WBDMIN       =  initial minimum accepted weight
C                                    for bond construction.
C                    WSTMAX       =  initial maximum accepted weight
C                                    for antibond construction.
C                    WBDCRT       =  critical weight limit for bond
C                                    construction.
C                    WSTEP        =  weight stepping size to decrease/
C                                    increase the bond/antibond weight
C                                    limits.
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
C                    IVEC,XVEC    =  int/flp scratch array of vector
C                                    type.
C                    XMAT         =  flp scratch array of matrix type.
C                    BDATOM       =  NAO atomic index map.
C                    W            =  NAO weight vector in atomic order.
C                    C            =  NAO coefficient matrix in AO basis
C                                    with columns in atomic order.
C
C
C                  Output:
C
C                    NBOND        =  total # of bonds found in the NHB
C                                    part (that is, excluding the Core
C                                    and Rydberg 'bonds').
C                    WBOND (x)    =  will contain the weight acceptance
C                                    criterion for x-center bond
C                                    formation during the NHO search
C                                    process. The final values when
C                                    exiting this routine will be the
C                                    lowest for each x-center bond
C                                    that lead to the final complete
C                                    set of NHO.
C                    WSTAR (x)    =  will contain the weight acceptance
C                                    criterion for x-center antibond
C                                    formation during the NHO search
C                                    process. The final values when
C                                    exiting this routine will be the
C                                    highest for each x-center antibond
C                                    that lead to the final complete
C                                    set of NHO.
C                    BDNCEN (J)   =  # of atomic centers for J-th bond.
C                                    The bond order is NHB-bonds/NCB/NRB.
C                                    Note that the # of NHB-bonds might
C                                    be < NHB size.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDNBAS (J)   =  # of basis functions (NHOs) for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDOCC (J)    =  # of occupied levels for J-th bond.
C                                    The bond order is NHB-bonds/NCB/NRB.
C                                    Note that the # of NHB-bonds might
C                                    be < NHB size.
C                    BDATOM       =  NAO atomic index map in NHB/NCB/NRB
C                                    order.
C                    ANGNHO (J,L) =  NBAS x (MXSHELL+1) matrix
C                                    containing the sum of the square
C                                    of the NAO coefficients per L-type
C                                    angular momentum for the J-th NHO.
C                                    The NHO order is NHB/NCB/NRB.
C                    H (I,J)      =  MXNBA x NBAS matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB/NCB/NRB
C                                    NHOs in that order, with all NHOs
C                                    belonging to a specific atomic
C                                    center being grouped together.
C                    W            =  NAO weight vector in NHB/NCB/NRB
C                                    order. The NAO weights in the NHB
C                                    part will have been overwritten
C                                    by the pre-NHO weights to determine
C                                    the NHB NHOs.
C                    C            =  NAO coefficient matrix in AO basis
C                                    with columns in NHB/NCB/NRB order.
C
C
C                  !!! IMPORTANT CONCEPT ABOUT BONDS AS DEFINED HERE!!!
C
C                  The concept of a bond is here synonymous with atomic
C                  center connectivity. That means that although a
C                  specific atomic center connection might lead to many
C                  occupied 'bonds' in the usual chemists point of view,
C                  all those occupied 'bonds' are treated here as one
C                  bond and its different occupied levels. As an example
C                  consider benzene C6H6. If BONDSIZE is set to 6, it
C                  will result in one! pi-bond involving the six carbon
C                  centers and this bond will have a total of 6 levels
C                  of which 3 are ocupied.
C
C                  Note also that the Core and Rydberg NHOs found here
C                  will be considered as bonds belonging to one atom
C                  with respective occupation level numbers 1 and 0.
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

         CHARACTER*16  STAGE

         LOGICAL     FAILED
         LOGICAL     LTRG,UTRG
         LOGICAL     MORE
         LOGICAL     NAO2NHO
         LOGICAL     RYD2HYB
         LOGICAL     SAVEP

         INTEGER     ATOM
         INTEGER     BONDSIZE,NBOSIZE
         INTEGER     I,J
         INTEGER     INCSIZE
         INTEGER     LTYPE
         INTEGER     MXCHOOSE
         INTEGER     MXNHBA,MXNBA,MXSHELL,MXCOL
         INTEGER     NAO,NHO
         INTEGER     NATOM,NHATOM
         INTEGER     NBAS
         INTEGER     NBOND
         INTEGER     NHB,NHA,NCB,NCA,NRB,NRA
         INTEGER     NHBA,NVBA,NCBA,NRBA
         INTEGER     NHCEN
         INTEGER     NSTEP,STEP
         INTEGER     OFF

         INTEGER     ATNCB   (1:NCA   )
         INTEGER     ATNHB   (1:NHA   )
         INTEGER     ATNRB   (1:NRA   )
         INTEGER     ATCIDX  (1:NCA   )
         INTEGER     ATRIDX  (1:NRA   )
         INTEGER     ATCOFF  (1:NCA   )
         INTEGER     ATHOFF  (1:NHA   )
         INTEGER     ATROFF  (1:NRA   )
         INTEGER     ATHCEN  (1:NHA   )
         INTEGER     ATHVAL  (1:NHA   )
         INTEGER     ATORD   (1:NHA   )
         INTEGER     BDATOM  (1:NBAS  )
         INTEGER     BDNBAS  (1:NBAS  )
         INTEGER     BDNCEN  (1:NBAS  )
         INTEGER     BDOCC   (1:NBAS  )
         INTEGER     COLMAP  (1:NBAS  )
         INTEGER     HSHELL  (1:NHB   )
         INTEGER     CSHELL  (1:NCB   )
         INTEGER     RSHELL  (1:NRB   )
         INTEGER     IVEC    (1:2*NBAS+2*NHB+2*NHA)
         INTEGER     NCHOOSE (1:BONDSIZE)
         INTEGER     NHYB    (1:NHA   )
         INTEGER     NWSTEP  (1:BONDSIZE)

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS)
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS)

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE,1:BONDSIZE)

         DOUBLE PRECISION  CSQR
         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  WBDMIN,WSTMAX,WBDCRT,WSTEP
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  W     (1:NBAS)
         DOUBLE PRECISION  WBOND (1:BONDSIZE)
         DOUBLE PRECISION  WSTAR (1:BONDSIZE)
         DOUBLE PRECISION  XVEC  (1:3*NBAS)

         DOUBLE PRECISION  ANGNHO (1:NBAS  ,0:MXSHELL)
         DOUBLE PRECISION  BOMAT  (1:NATOM ,1:NATOM  )
         DOUBLE PRECISION  C      (1:NBAS  ,1:NBAS   )
         DOUBLE PRECISION  P      (1:NBAS  ,1:NBAS   )
         DOUBLE PRECISION  H      (1:MXNBA ,1:NBAS   )
         DOUBLE PRECISION  PH     (1:NHB   ,1:NHB    )
         DOUBLE PRECISION  PHDEP  (1:NHB   ,1:NHB    )
         DOUBLE PRECISION  PHSUB  (1:NHB   ,1:NHB    )
         DOUBLE PRECISION  SAH    (1:MXNHBA,1:MXNHBA )
         DOUBLE PRECISION  XMAT   (1:NBAS  ,1:MXCOL  )

         PARAMETER  (ZERO  = 0.D0 )
C
C
C------------------------------------------------------------------------
C
C
C             ...reorder the NAO weights, the atomic index map and the
C                coefficient matrix from atomic to NHB/NCB/NRB order.
C
C
C         CALL    MAT__PRINT_V_INTEGER_NOZEROS
C     +
C     +                ( 1,
C     +                  ' COLMAP ',
C     +                  NBAS,
C     +                  NBAS,
C     +                  COLMAP )
C     +
C     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' Hybrid atomic ordering vector ',
     +                  1,NHA,
     +                  1,NHA,
     +                  ATORD )
     +
     +
         CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +              ( NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,
     +                COLMAP,
     +                XVEC,
     +
     +                        W )
     +
     +
C         CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
C     +
C     +                ( 6,
C     +                  ' Atomic reordered W vector ',
C     +                  NBAS,
C     +                  NBAS,
C     +                  W )
C     +
C     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,
     +                COLMAP,
     +                IVEC,
     +
     +                        BDATOM )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,NBAS,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        C )
     +
     +
C
C
C             ...initialize the # of bond/antibond weight limit steps.
C
C
         DO NHCEN = 1,BONDSIZE
            NWSTEP (NHCEN) = 1
         END DO

         INCSIZE = 0
C
C
C             ...start iterates on NHO search. Increase # of
C                bond/antibond weight limit steps for appropriate
C                bond size by + 1.
C
C
 1000    IF (INCSIZE.GT.BONDSIZE) THEN
             INCSIZE = 1
             NWSTEP (1) = NWSTEP (1) + 1
         ELSE IF (INCSIZE.GT.0) THEN
             NWSTEP (INCSIZE) = NWSTEP (INCSIZE) + 1
         END IF

         STAGE = 'start'

         CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +              ( STAGE,
     +                BONDSIZE,
     +                MAXOCC,
     +                NHCEN,
     +                NWSTEP,
     +                WBDMIN,WSTMAX,WSTEP,
     +                WBOND,WSTAR )
     +
     +
         LTRG  = .TRUE.
         UTRG  = .FALSE.
         SAVEP = .TRUE.

         CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +              ( NBAS,NHB,
     +                NBAS,NBAS,
     +                NBAS,NHB,
     +                NHB,
     +                NHB,NBAS,
     +                0,
     +                SAVEP,LTRG,UTRG,
     +                C,
     +                P,
     +                XVEC,
     +
     +                        XMAT )
     +
     +
         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NHB,
     +                NHB,NHB,
     +                NHB,NHB,
     +                XMAT,
     +
     +                        PH )
     +
     +
C         CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
C     +
C     +              ( NBAS,NHB,NHB,
C     +                NHB,
C     +                .TRUE.,
C     +
C     +                        XVEC,
C     +                        XMAT )
C     +
C     +
C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' PH matrix ',
C     +                  NHB,NHB,
C     +                  NHB,0,
C     +                  PH )
C     +
C     +
C         CALL    MAT__PRINT_V_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' Eigenvalues of PH matrix ',
C     +                  NHB,
C     +                  NHB,
C     +                  XVEC )
C     +
C     +
C
C
C             ...initialize search for pre-NHOs.
C
C
         CALL  MAT__W_EQ_ZERO_INTEGER  (NBAS,NBAS,BDNBAS)
         CALL  MAT__W_EQ_ZERO_INTEGER  (NBAS,NBAS,BDNCEN)
         CALL  MAT__W_EQ_ZERO_INTEGER  (NBAS,NBAS,BDOCC)
         CALL  MAT__C_EQ_ZERO_INTEGER  (NBOSIZE,NBAS,NBOSIZE,NBAS,BDBAS)
         CALL  MAT__C_EQ_ZERO_INTEGER  (NBOSIZE,NBAS,NBOSIZE,NBAS,BDCEN)
         CALL  MAT__C_EQ_ZERO_FLOAT    (MXNBA,NBAS,MXNBA,NBAS,H)
C
C
C             ...perform the loop over all allowed atomic hybrid
C                center combinations. Start with all hybrid atoms
C                in their optimum order using the optimum ordering
C                array. Copy initial PH matrix to the matrix that
C                will accumulate all depletions during a specific
C                bondsize run (if any). Keep going in reducing
C                the weight limits for each x-centered bond type
C                and check for pre-NHOs as long as that is still
C                possible.
C
C
         NBOND = 0

         NHATOM = 0
         DO I = 1,NHA
            ATOM = ATORD (I)
            NHATOM = NHATOM + 1
            ATHCEN (NHATOM) = ATOM
            NHYB (ATOM) = 0
         END DO

         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NHB,NHB,
     +                NHB,NHB,
     +                NHB,NHB,
     +                PH,
     +
     +                        PHDEP )
     +
     +
         DO NHCEN = 1,BONDSIZE

            NSTEP = NWSTEP (NHCEN)
            WBOND (NHCEN) = WBDMIN
            WSTAR (NHCEN) = WSTMAX

            DO STEP = 1,NSTEP

               IF (NHATOM.GT.0 .AND. NHCEN.LE.NHATOM) THEN

                   CALL  NLO__CHECK_NHO_ALL_X_CENTERS
     +
     +                        ( NATOM,
     +                          BONDSIZE,NBOSIZE,
     +                          NHCEN,NHATOM,
     +                          MXNHBA,MXNBA,
     +                          ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                          RYD2HYB,
     +                          IVEC (1),IVEC (NHA+1),
     +                          NHB,NHA,
     +                          WBOND (NHCEN),WSTAR (NHCEN),
     +                          SAH,
     +                          PH,PHSUB,
     +                          BOMAT,
     +                          MXCHOOSE,NCHOOSE (NHCEN),
     +                          CHOOSE (1,1,NHCEN),
     +                          IVEC (2*NHA+1),
     +                          XVEC,XMAT,
     +
     +                                   FAILED,
     +                                   NBOND,
     +                                   NHYB,
     +                                   BDNCEN,
     +                                   BDCEN,
     +                                   BDNBAS,
     +                                   BDBAS,
     +                                   BDOCC,
     +                                   W,H,
     +                                   PHDEP )
     +
     +
                   IF (FAILED) THEN
                       STAGE = 'failure'
                       CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                       ( STAGE,
     +                         BONDSIZE,
     +                         MAXOCC,
     +                         NHCEN,
     +                         NWSTEP,
     +                         WBDMIN,WSTMAX,WSTEP,
     +                         WBOND,WSTAR )
     +
     +
                       INCSIZE = INCSIZE + 1
                       GOTO 1000
                   END IF
C
C
C             ...check, if any hybrid atoms have a nonexhaustive bond
C                forming pre-NHO set. Determine their number and check
C                for more bond forming pre-NHO's. If there are more
C                allowed, update the PH matrix with all depletions
C                accumulated during the present bondsize run and
C                continue.
C
C
                   NHATOM = 0
                   DO I = 1,NHA
                      ATOM = ATORD (I)
                      IF (NHYB (ATOM) .LT. ATHVAL (ATOM)) THEN
                          NHATOM = NHATOM + 1
                          ATHCEN (NHATOM) = ATOM
                      END IF
                   END DO

                   CALL  MAT__C_EQ_A_FLOAT
     +
     +                        ( NHB,NHB,
     +                          NHB,NHB,
     +                          NHB,NHB,
     +                          PHDEP,
     +
     +                               PH )
     +
     +
               END IF

               IF (STEP.NE.NSTEP) THEN
                   WBOND (NHCEN) = WBOND (NHCEN) - WSTEP
                   WSTAR (NHCEN) = WSTAR (NHCEN) + WSTEP
               END IF

            END DO

         END DO
C
C
C             ...try to complete the present sets of bond forming
C                pre-NHOs on each atom with eventual Empty-pair
C                pre-NHOs and possibly Rydberg pre-NHOs. If the
C                routine is unable to do so, we have to retry the
C                the bond forming pre-NHO search with lower bond
C                weight limits. If the routine is successful, the
C                complete set of orthogonalized NHOs on each atom
C                have been obtained.
C
C
 2000    STAGE = 'complete'

         CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +              ( STAGE,
     +                BONDSIZE,
     +                MAXOCC,
     +                NHCEN,
     +                NWSTEP,
     +                WBDMIN,WSTMAX,WSTEP,
     +                WBOND,WSTAR )
     +
     +
         CALL  NLO__COMPLETE_NHO_SPACE
     +
     +              ( NBAS,NHATOM,
     +                NBOSIZE,
     +                MXNHBA,MXNBA,
     +                ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                IVEC (1),
     +                NHB,NHA,
     +                NHYB,
     +                WBDCRT,WSTAR (1),
     +                SAH,
     +                P,PH,PHSUB,
     +                W,C,
     +                IVEC (NHATOM+1),
     +                XVEC,XMAT,
     +
     +                        FAILED,
     +                        MORE,
     +                        NBOND,
     +                        BDNCEN,
     +                        BDCEN,
     +                        BDNBAS,
     +                        BDBAS,
     +                        BDOCC,
     +                        H )
     +
     +
         IF (FAILED) THEN
             STAGE = 'failure complete'
             CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                  ( STAGE,
     +                    BONDSIZE,
     +                    MAXOCC,
     +                    NHCEN,
     +                    NWSTEP,
     +                    WBDMIN,WSTMAX,WSTEP,
     +                    WBOND,WSTAR )
     +
     +
             INCSIZE = INCSIZE + 1
             GOTO 1000
         END IF

         IF (MORE) THEN

             STAGE = 'more'

             CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                  ( STAGE,
     +                    BONDSIZE,
     +                    MAXOCC,
     +                    NHCEN,
     +                    NWSTEP,
     +                    WBDMIN,WSTMAX,WSTEP,
     +                    WBOND,WSTAR )
     +
     +
             DO NHCEN = 1,BONDSIZE

                WBOND (NHCEN) = WBOND (NHCEN) - WSTEP
                WSTAR (NHCEN) = WSTAR (NHCEN) + WSTEP

                IF (NHATOM.GT.0 .AND. NHCEN.LE.NHATOM) THEN

                    CALL  NLO__CHECK_NHO_ALL_X_CENTERS
     +
     +                         ( NATOM,
     +                           BONDSIZE,NBOSIZE,
     +                           NHCEN,NHATOM,
     +                           MXNHBA,MXNBA,
     +                           ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                           RYD2HYB,
     +                           IVEC (1),IVEC (NHA+1),
     +                           NHB,NHA,
     +                           WBOND (NHCEN),WSTAR (NHCEN),
     +                           SAH,
     +                           PH,PHSUB,
     +                           BOMAT,
     +                           MXCHOOSE,NCHOOSE (NHCEN),
     +                           CHOOSE (1,1,NHCEN),
     +                           IVEC (2*NHA+1),
     +                           XVEC,XMAT,
     +
     +                                    FAILED,
     +                                    NBOND,
     +                                    NHYB,
     +                                    BDNCEN,
     +                                    BDCEN,
     +                                    BDNBAS,
     +                                    BDBAS,
     +                                    BDOCC,
     +                                    W,H,
     +                                    PHDEP )
     +
     +
                    IF (FAILED) THEN
                        STAGE = 'failure'
                        CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                        ( STAGE,
     +                          BONDSIZE,
     +                          MAXOCC,
     +                          NHCEN,
     +                          NWSTEP,
     +                          WBDMIN,WSTMAX,WSTEP,
     +                          WBOND,WSTAR )
     +
     +
                        INCSIZE = INCSIZE + 1
                        GOTO 1000
                    END IF

                    NHATOM = 0
                    DO I = 1,NHA
                       ATOM = ATORD (I)
                       IF (NHYB (ATOM) .LT. ATHVAL (ATOM)) THEN
                           NHATOM = NHATOM + 1
                           ATHCEN (NHATOM) = ATOM
                       END IF
                    END DO

                    CALL  MAT__C_EQ_A_FLOAT
     +
     +                         ( NHB,NHB,
     +                           NHB,NHB,
     +                           NHB,NHB,
     +                           PHDEP,
     +
     +                                PH )
     +
     +
                END IF
            END DO
            GOTO 2000
         END IF

         STAGE = 'success'

         CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +              ( STAGE,
     +                BONDSIZE,
     +                MAXOCC,
     +                NHCEN,
     +                NWSTEP,
     +                WBDMIN,WSTMAX,WSTEP,
     +                WBOND,WSTAR )
     +
     +
C
C
C             ...form now the Core and Rydberg NHO orbitals.
C
C
         NAO2NHO = .TRUE.

         CALL  NLO__FORM_ATOMIC_NONVALENCE_NHO
     +
     +              ( NBAS,MXNBA,
     +                NBOSIZE,
     +                ATNCB,ATCIDX,ATCOFF,
     +                NCB,NCA,
     +                1,NHB,
     +                NAO2NHO,
     +                P,
     +                C (1,NHB+1),
     +                IVEC (1),IVEC (MXNBA+1),
     +                XVEC,
     +                XMAT,
     +
     +                        BDNCEN (NHB+1),
     +                        BDCEN (1,NHB+1),
     +                        BDNBAS (NHB+1),
     +                        BDBAS (1,NHB+1),
     +                        BDOCC (NHB+1),
     +                        H (1,NHB+1))
     +
     +
         CALL  NLO__FORM_ATOMIC_NONVALENCE_NHO
     +
     +              ( NBAS,MXNBA,
     +                NBOSIZE,
     +                ATNRB,ATRIDX,ATROFF,
     +                NRB,NRA,
     +                0,NHB+NCB,
     +                NAO2NHO,
     +                P,
     +                C (1,NHB+NCB+1),
     +                IVEC (1),IVEC (MXNBA+1),
     +                XVEC,
     +                XMAT,
     +
     +                        BDNCEN (NHB+NCB+1),
     +                        BDCEN (1,NHB+NCB+1),
     +                        BDNBAS (NHB+NCB+1),
     +                        BDBAS (1,NHB+NCB+1),
     +                        BDOCC (NHB+NCB+1),
     +                        H (1,NHB+NCB+1))
     +
     +
         CALL    MAT__PRINT_V_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' BDNBAS vector ',
     +                  NBAS,
     +                  NBAS,
     +                  BDNBAS )
     +
     +
         CALL    MAT__PRINT_V_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' BDNCEN vector ',
     +                  NBAS,
     +                  NBAS,
     +                  BDNCEN )
     +
     +
         CALL    MAT__PRINT_V_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' BDOCC vector ',
     +                  NBAS,
     +                  NBAS,
     +                  BDOCC )
     +
     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' BDCEN array ',
     +                  NBOSIZE,NBAS,
     +                  NBOSIZE,NBAS,
     +                  BDCEN )
     +
     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' BDBAS array ',
     +                  NBOSIZE,NBAS,
     +                  NBOSIZE,NBAS,
     +                  BDBAS )
     +
     +
C
C
C             ...calculate the matrix containing the sum of the
C                square of the NAO coefficients per L-type angular
C                momentum for all NHOs. The NHOs are analyzed in
C                the following sequence; valence, core, hybrid
C                Rydberg and Rydberg.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NBAS,MXSHELL+1,
     +                NBAS,MXSHELL+1,
     +
     +                        ANGNHO )
     +
     +
         DO ATOM = 1,NHA
            OFF = ATHOFF (ATOM)
            NHBA = ATNHB (ATOM)
            NVBA = ATHVAL (ATOM)
            DO I = 1,NHBA
               NAO = OFF + I
               LTYPE = HSHELL (NAO)
               DO J = 1,NVBA
                  NHO = OFF + J
                  CSQR = H (I,NHO) ** 2
                  ANGNHO (NHO,LTYPE) = ANGNHO (NHO,LTYPE) + CSQR
               END DO
            END DO
         END DO

         DO ATOM = 1,NCA
            OFF = ATCOFF (ATOM)
            NCBA = ATNCB (ATOM)
            DO I = 1,NCBA
               NAO = OFF + I
               LTYPE = CSHELL (NAO)
               DO J = 1,NCBA
                  NHO = NHB + OFF + J
                  CSQR = H (I,NHO) ** 2
                  ANGNHO (NHO,LTYPE) = ANGNHO (NHO,LTYPE) + CSQR
               END DO
            END DO
         END DO

         IF (RYD2HYB) THEN
             DO ATOM = 1,NHA
                OFF = ATHOFF (ATOM)
                NHBA = ATNHB (ATOM)
                NVBA = ATHVAL (ATOM)
                NRBA = NHBA - NVBA
                DO I = 1,NHBA
                   NAO = OFF + I
                   LTYPE = HSHELL (NAO)
                   DO J = 1,NRBA
                      NHO = OFF + NVBA + J
                      CSQR = H (I,NHO) ** 2
                      ANGNHO (NHO,LTYPE) = ANGNHO (NHO,LTYPE) + CSQR
                   END DO
                END DO
             END DO
         END IF

         DO ATOM = 1,NRA
            OFF = ATROFF (ATOM)
            NRBA = ATNRB (ATOM)
            DO I = 1,NRBA
               NAO = OFF + I
               LTYPE = RSHELL (NAO)
               DO J = 1,NRBA
                  NHO = NHB + NCB + OFF + J
                  CSQR = H (I,NHO) ** 2
                  ANGNHO (NHO,LTYPE) = ANGNHO (NHO,LTYPE) + CSQR
               END DO
            END DO
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
