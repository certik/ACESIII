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
         SUBROUTINE  NLO__GENER_NBO_ORBITALS
     +
     +                    ( NBAS,NATOM,NBOND,
     +                      BONDSIZE,NBOSIZE,
     +                      MXNBA,MXSHELL,
     +                      NHB,NCB,NRB,
     +                      NHA,NCA,NRA,
     +                      ATNHB,ATHOFF,
     +                      ATNCB,ATCIDX,ATCOFF,
     +                      ATNRB,ATRIDX,ATROFF,
     +                      BASBEG,BASEND,
     +                      COLMAP,
     +                      WBOND,WSTAR,
     +                      WRYD,WOCC,
     +                      BDNCEN,BDCEN,
     +                      BDNBAS,BDBAS,
     +                      BDOCC,
     +                      ANGNHO,
     +                      P,H,
     +                      PH,
     +                      IVEC,
     +                      XVEC,XMAT,
     +
     +                              NLA,NEA,NYA,
     +                              NLB,NBB,NEB,NAB,NYB,
     +                              ATNLB,ATLIDX,ATLOFF,
     +                              ATNEB,ATEIDX,ATEOFF,
     +                              ATNYB,ATYIDX,ATYOFF,
     +                              NBOBD,
     +                              ANGNBO,
     +                              BDSIZE,BDATOM,
     +                              B,W,Q,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NBO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the sets of natural bond
C                orbitals (NBOs) by combining the final set of NHOs
C                in the appropriate way.
C
C
C                  Input:
C
C                    NBAS         =  total # of AOs in AO basis
C                    NATOM        =  total # of atomic centers
C                    NBOND        =  total # of bonds found in the NHB
C                                    part.
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    MXSHELL      =  largest l-shell value
C                    NxB          =  total # of Hybrid, Core and
C                                    Rydberg NAOs (x=H,C,R)
C                    NxA          =  total # of hybrid, core and
C                                    Rydberg atoms (x=H,C,R)
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid
C                                    NAOs on all hybrid atoms preceeding
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
C                    BASBEG (A)   =  first basis index number for atom A
C                    BASEND (A)   =  last basis index number for atom A
C                    COLMAP (I)   =  will contain the column map in the
C                                    active sense for the pure NBO
C                                    ordering NHB -> NLB/NBB/NEB/NAB/NYB
C                                    and the subsequent mixed NBO/NAO
C                                    ordering NHB/NCB -> NCB/NHB to bring
C                                    the Core NAOs up front. COLMAP (I)
C                                    contains the position index of the
C                                    I-th NBO in the final respective
C                                    order.
C                    WBOND (x)    =  lowest weight criterion for each
C                                    x-center bond.
C                    WSTAR (x)    =  highest weight criterion for each
C                                    x-center antibond.
C                    WRYD         =  highest weight below which a one-
C                                    center bond is classified as an
C                                    atomic Rydberg NBO.
C                    WOCC         =  the weight limit to decide when
C                                    a NBO is considered occupied or
C                                    virtual.
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
C                    ANGNHO (J,L) =  NBAS x (MXSHELL+1) matrix
C                                    containing the sum of the square
C                                    of the NAO coefficients per L-type
C                                    angular momentum for the J-th NHO.
C                                    The NHO order is NHB/NCB/NRB.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    H (I,J)      =  MXNBA x NBAS matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic Hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHOs in
C                                    NHB/NCB/NRB order, with all NHOs
C                                    in each section belonging to a
C                                    specific atomic center being
C                                    grouped together.
C                    PH           =  will contain the NHB x NHB valence
C                                    occupation matrix in NAO basis.
C                    IVEC         =  int scratch array of vector type.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    BDSIZE       =  full NAO bond sizes (all = 1).
C                    BDATOM       =  full NAO atomic index map with
C                                    1st column in NHB/NCB/NRB order.
C                    W            =  weight vector in NHB/NCB/NRB order.
C                                    The NCB and NRB parts contain the
C                                    original Core and Rydberg NAO
C                                    weights. The NHB part was used
C                                    during NHO generation and will
C                                    be replaced by the NBO Lone-pair,
C                                    Bond, Empty-pair, Antibond and
C                                    Rydberg weights.
C                    C            =  full NAO coefficient matrix in AO
C                                    basis with columns in NHB/NCB/NRB
C                                    order.
C
C
C                  Output:
C
C                    NxA          =  total # of Lone-pair, Empty-pair
C                                    and Rydberg atoms found (x=L,E,Y).
C                    NxB          =  total # of Lone-pair, Bond,
C                                    Empty-pair, Antibond and Rydberg
C                                    NBOs found (x=L,B,E,A,Y).
C                    ATNLB (A)    =  # of Lone-pair NBOs on Lone-pair
C                                    atom A.
C                    ATLIDX (A)   =  atomic index for Lone-pair atom A.
C                    ATLOFF (A)   =  index offset for Lone-pair NBOs for
C                                    Lone-pair atom A. This index is
C                                    equal to the total number of
C                                    Lone-pair NBOs on all Lone-pair
C                                    atoms preceeding Lone-pair atom A.
C                    ATNEB (A)    =  # of Empty-pair NBOs on Empty-pair
C                                    atom A.
C                    ATEIDX (A)   =  atomic index for Empty-pair atom A.
C                    ATEOFF (A)   =  index offset for Empty-pair NBOs for
C                                    Empty-pair atom A. This index is
C                                    equal to the total number of
C                                    Empty-pair NBOs on all Empty-pair
C                                    atoms preceeding Empty-pair atom A.
C                    ATNYB (A)    =  # of Rydberg NBOs on Rydberg atom A.
C                    ATYIDX (A)   =  atomic index for Rydberg atom A.
C                    ATYOFF (A)   =  index offset for Rydberg NBOs for
C                                    Rydberg atom A. This index is
C                                    equal to the total number of
C                                    Rydberg NBOs on all Rydberg atoms
C                                    preceeding Rydberg atom A.
C                    NBOBD (I)    =  contains the bond index number for
C                                    the I-th NBO. This array is the
C                                    handle for accessing and modifying
C                                    info of the bonds sitting in the
C                                    arrays BDNCEN,BDCEN,BDBAS and
C                                    BDOCC. The NBOBD elements are in
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB order.
C                    ANGNBO (J,L) =  NBAS x (MXSHELL+1) matrix containing
C                                    the sum of the square of the NAO
C                                    coefficients per L-type angular
C                                    momentum for the J-th NBO. The NBO
C                                    order is NCB/NLB/NBB/NEB/NAB/NYB/
C                                    NRB.
C                    BDSIZE       =  full NBO bond sizes in NBO order
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB.
C                    BDATOM       =  full NBO atomic index map in NBO
C                                    order NCB/NLB/NBB/NEB/NAB/NYB/NRB.
C                    B (I,J)      =  NBOSIZE x NBAS matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO. The NBO
C                                    column index is in NCB/NLB/NBB/NEB/
C                                    NAB/NYB/NRB order.
C                    W            =  NBO weight vector with elements in
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB order.
C                    Q            =  NBO occupation interaction order
C                                    vector with elements in NCB/NLB/
C                                    NBB/NEB/NAB/NYB/NRB order. These
C                                    interaction orders are given by
C                                    the square root of the sum of the
C                                    square offdiagonal occupation
C                                    matrix elements:
C
C                                      Q (I) = sqrt ( sum   P (I,K)**2 )
C                                                    K.ne.I
C
C                                    They will be used here to determine
C                                    degenerate NBOs during the NBO
C                                    symmetrization procedure.
C                    C            =  full NBO coefficient matrix
C                                    in AO basis with columns in
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB order.
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

         LOGICAL     LTRG,UTRG
         LOGICAL     SAVEP

         INTEGER     ATOM
         INTEGER     ATNEW,ATOLD
         INTEGER     BOND
         INTEGER     BONDSIZE,NBOSIZE
         INTEGER     I,J
         INTEGER     INDEX
         INTEGER     LTYPE
         INTEGER     MXNBA,MXSHELL
         INTEGER     NBAS,NATOM
         INTEGER     NBO,NBOOLD
         INTEGER     NBOND
         INTEGER     NHB,NHA,NCB,NCA,NRB,NRA
         INTEGER     NHCEN
         INTEGER     NHO,NNHO
         INTEGER     NLB,NLA,NBB,NEB,NEA,NAB,NYB,NYA
         INTEGER     NLBA,NEBA,NYBA
         INTEGER     NLBSUM,NEBSUM,NYBSUM
         INTEGER     NMOVE
         INTEGER     NOCC
         INTEGER     OFFLP,OFFBD,OFFEP,OFFAB,OFFRY

         INTEGER     ATNCB   (1:NCA      )
         INTEGER     ATNHB   (1:NHA      )
         INTEGER     ATNRB   (1:NRA      )
         INTEGER     ATNLB   (1:NATOM    )
         INTEGER     ATNEB   (1:NATOM    )
         INTEGER     ATNYB   (1:NATOM    )
         INTEGER     ATCIDX  (1:NCA      )
         INTEGER     ATRIDX  (1:NRA      )
         INTEGER     ATLIDX  (1:NATOM    )
         INTEGER     ATEIDX  (1:NATOM    )
         INTEGER     ATYIDX  (1:NATOM    )
         INTEGER     ATCOFF  (1:NCA      )
         INTEGER     ATHOFF  (1:NHA      )
         INTEGER     ATROFF  (1:NRA      )
         INTEGER     ATLOFF  (1:NATOM    )
         INTEGER     ATEOFF  (1:NATOM    )
         INTEGER     ATYOFF  (1:NATOM    )
         INTEGER     BASBEG  (1:NATOM    )
         INTEGER     BASEND  (1:NATOM    )
         INTEGER     BDNBAS  (1:NBAS     )
         INTEGER     BDNCEN  (1:NBAS     )
         INTEGER     BDOCC   (1:NBAS     )
         INTEGER     BDSIZE  (1:NBAS     )
         INTEGER     COLMAP  (1:NCB+NHB  )
         INTEGER     IVEC    (1:NBAS     )
         INTEGER     NBOBD   (1:NBAS     )

         INTEGER     BDATOM  (1:NBAS   ,1:BONDSIZE)
         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS    )
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS    )

         DOUBLE PRECISION  ANGVAL
         DOUBLE PRECISION  CSQR
         DOUBLE PRECISION  WRYD,WOCC
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  Q     (1:NBAS     )
         DOUBLE PRECISION  W     (1:NBAS     )
         DOUBLE PRECISION  WBOND (1:BONDSIZE )
         DOUBLE PRECISION  WSTAR (1:BONDSIZE )
         DOUBLE PRECISION  XVEC  (1:NBAS+NBAS)

         DOUBLE PRECISION  ANGNBO (1:NBAS   ,0:MXSHELL)
         DOUBLE PRECISION  ANGNHO (1:NBAS   ,0:MXSHELL)
         DOUBLE PRECISION  B      (1:NBOSIZE,1:NBAS   )
         DOUBLE PRECISION  C      (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  P      (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  H      (1:MXNBA  ,1:NBAS   )
         DOUBLE PRECISION  PH     (1:NHB    ,1:NHB    )
         DOUBLE PRECISION  XMAT   (1:NBAS   ,1:NBAS   )

         PARAMETER  (ZERO  = 0.D0 )
         PARAMETER  (ONE   = 1.D0 )
C
C
C------------------------------------------------------------------------
C
C
C             ...form the full NHB x NHB occupation matrix PH.
C
C
         LTRG  = .TRUE.
         UTRG  = .TRUE.
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
C
C
C             ...loop over all NHO bonds found and accumulate the NBO
C                weights, the coefficient matrix, the NBO bond sizes,
C                the NBO atomic index map, the NHO bond index
C                numbers and the local column NBO ordering map
C                NHB -> NLB/NBB/NEB/NAB/NYB with added offsets for
C                temporary protection. The difference between the
C                offsets is set to NHB, since this ensures no overlap
C                of indices between different types of NBOs.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NBOSIZE,NBAS,
     +                NBOSIZE,NBAS,
     +
     +                        B )
     +
     +
         NBO = 0
         NLB = 0
         NBB = 0
         NEB = 0
         NAB = 0
         NYB = 0
         NBOOLD = 0

         OFFLP = 0
         OFFBD = OFFLP + NHB
         OFFEP = OFFBD + NHB
         OFFAB = OFFEP + NHB
         OFFRY = OFFAB + NHB

         DO BOND = 1,NBOND

            NOCC = BDOCC (BOND)
            NHCEN = BDNCEN (BOND)

            CALL  NLO__FORM_NBO_X_CENTERS
     +
     +              ( NBAS,NHCEN,NOCC,
     +                BONDSIZE,NBOSIZE,
     +                MXNBA,
     +                ATNHB,ATHOFF,
     +                NHB,NHA,
     +                OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,
     +                BDCEN (1,BOND),
     +                BDBAS (1,BOND),
     +                WBOND (NHCEN),
     +                WSTAR (NHCEN),
     +                WRYD,WOCC,
     +                C,
     +                PH,
     +                H,
     +                IVEC,
     +
     +                        NBO,
     +                        NLB,NBB,NEB,NAB,NYB,
     +                        COLMAP,
     +                        B (1,NBOOLD+1),
     +                        W (NBOOLD+1),
     +                        BDSIZE,BDATOM,
     +                        XMAT (1,NBO+1) )
     +
     +
            DO I = NBOOLD+1,NBO
               NBOBD (I) = BOND
            END DO

            NBOOLD = NBO
         END DO
C
C
C             ...check, if the # of Lone-pairs, Bonds, Empty-pairs,
C                Antibonds and Rydbergs add up ok.
C
C
         IF ((NLB+NBB+NEB+NAB+NYB).NE.NHB) THEN
             WRITE (*,*) ' Problems in finding correct # of NBOs! '
             WRITE (*,*) ' NLB,NBB,NEB,NAB,NYB,NHB = ',
     +                     NLB,NBB,NEB,NAB,NYB,NHB
             WRITE (*,*) ' nlo__gener_nbo_orbitals '
             WRITE (1,*) ' Problems in finding correct # of NBOs! '
             WRITE (1,*) ' NLB,NBB,NEB,NAB,NYB,NHB = ',
     +                     NLB,NBB,NEB,NAB,NYB,NHB
             WRITE (1,*) ' nlo__gener_nbo_orbitals '
             STOP
         END IF
C
C
C             ...complete the NHO bond index number array by adding
C                the Core and Rydberg NHO bonds. This is just a one
C                to one index mapping. Complete also the NBO expansion
C                coefficient matrix, which simply contains a +1 in its
C                first row for the Core and Rydberg NBOs, as these
C                consist of just one NHO at this stage.
C
C
         DO I = NHB+1,NBAS
            NBOBD (I) = I
            B (1,I) = ONE
         END DO
C
C
C             ...copy the NHB section of the NBO coefficient
C                accumulation matrix to the final NBO coefficient
C                matrix.
C
C
         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NHB,
     +                NBAS,NBAS,
     +                NBAS,NHB,
     +                XMAT,
     +
     +                        C )
     +
     +
C
C
C             ...remove the protective offsets on the local column
C                NBO ordering NHB -> NLB/NBB/NEB/NAB/NYB and reorder
C                all the NHB related data:
C
C                    i) the NBO weights
C                   ii) the NBO coefficients in AO basis
C                  iii) the NBO -> NHO expansion coefficient matrix
C                   iv) the vector containing the NHO bond index numbers
C                    v) the NBO bond size vector
C                   vi) the NBO atomic index array
C
C                Use the obsolete PH matrix for scratch.
C
C
         DO I = 1,NHB
            INDEX = COLMAP (I)
            IF (INDEX.GT.OFFRY) THEN
                INDEX = INDEX - OFFRY + NLB + NBB + NEB + NAB
            ELSE IF (INDEX.GT.OFFAB) THEN
                INDEX = INDEX - OFFAB + NLB + NBB + NEB
            ELSE IF (INDEX.GT.OFFEP) THEN
                INDEX = INDEX - OFFEP + NLB + NBB
            ELSE IF (INDEX.GT.OFFBD) THEN
                INDEX = INDEX - OFFBD + NLB
            ELSE IF (INDEX.GT.OFFLP) THEN
                INDEX = INDEX - OFFLP
            END IF
            COLMAP (I) = INDEX
         END DO

         CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +              ( NHB,
     +                NHB,
     +                NHB,
     +                NHB,
     +                COLMAP,
     +                XVEC,
     +
     +                        W )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NHB,
     +                NHB,
     +                NHB,
     +                NBAS,
     +                NBAS,NHB,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        C )
     +
     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NHB,
     +                NHB,
     +                NHB,
     +                NHB,
     +                COLMAP,
     +                IVEC,
     +
     +                        NBOBD )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBOSIZE,NHB,
     +                NHB,
     +                NHB,
     +                NBOSIZE,
     +                NBOSIZE,NHB,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        B )
     +
     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NHB,
     +                NHB,
     +                NHB,
     +                NHB,
     +                COLMAP,
     +                IVEC,
     +
     +                        BDSIZE )
     +
     +
         DO I = 1,BONDSIZE

            CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +                 ( NHB,
     +                   NHB,
     +                   NHB,
     +                   NHB,
     +                   COLMAP,
     +                   IVEC,
     +
     +                           BDATOM (1,I) )
     +
     +
         END DO
C
C
C             ...at this stage all NLB lone pairs are conglomerated,
C                however, due to their finding process, they are
C                not ordered according to the atomic index. Lone
C                pairs corresponding to the same atom are scattered
C                through the set of NLB lone pairs. We thus generate
C                first the COLMAP array to bring the lone pairs into
C                atomic order and then reorder them together with
C                all the associated data. After that we analyze the
C                lone pair space and get some info from it to be used
C                for symmetrization.
C
C
         DO I = 1,NATOM
            IVEC (I) = 0
         END DO
         DO I = 1,NLB
            ATOM = BDATOM (I,1)
            OFFLP = ATOM * NLB
            IVEC (ATOM) = IVEC (ATOM) + 1
            COLMAP (I) = IVEC (ATOM) + OFFLP
         END DO
         DO I = NATOM,1,-1
            IVEC (I) = 0
            DO J = 1,I-1
               IVEC (I) = IVEC (I) + IVEC (J)
            END DO
         END DO
         DO I = 1,NLB
            INDEX = COLMAP (I)
            ATOM = INDEX / NLB
            OFFLP = ATOM * NLB
            COLMAP (I) = INDEX - OFFLP + IVEC (ATOM)
         END DO

         CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +              ( NLB,
     +                NLB,
     +                NLB,
     +                NLB,
     +                COLMAP,
     +                XVEC,
     +
     +                        W )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NLB,
     +                NLB,
     +                NLB,
     +                NBAS,
     +                NBAS,NLB,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        C )
     +
     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NLB,
     +                NLB,
     +                NLB,
     +                NLB,
     +                COLMAP,
     +                IVEC,
     +
     +                        NBOBD )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBOSIZE,NLB,
     +                NLB,
     +                NLB,
     +                NBOSIZE,
     +                NBOSIZE,NLB,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        B )
     +
     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NLB,
     +                NLB,
     +                NLB,
     +                NLB,
     +                COLMAP,
     +                IVEC,
     +
     +                        BDSIZE )
     +
     +
         DO I = 1,BONDSIZE

            CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +                 ( NLB,
     +                   NLB,
     +                   NLB,
     +                   NLB,
     +                   COLMAP,
     +                   IVEC,
     +
     +                           BDATOM (1,I) )
     +
     +
         END DO
C
C
C             ...all data related to the NHB section is now locally
C                ordered in the form NHB -> NLB/NBB/NEB/NAB/NYB.
C                Determine global mixed column NBO/NAO ordering NHB/NCB
C                -> NCB/NHB. Reorder the NBO coefficient matrix from
C                NHB/NCB order to NCB/NHB order, such that the Core
C                NBOs come first. Reorder also the NHB/NCB part of the
C                weight vector and the NHB/NCB parts of the NBO bond
C                size vector and the NBO atomic index array.
C
C
         DO I = 1,NHB
            COLMAP (I) = I + NCB
         END DO

         DO I = 1,NCB
            COLMAP (NHB+I) = I
         END DO

         NMOVE = NCB + NHB

         CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +              ( NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                COLMAP,
     +                XVEC,
     +
     +                        W )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NBAS,
     +                NBAS,NMOVE,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        C )
     +
     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                COLMAP,
     +                IVEC,
     +
     +                        NBOBD )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBOSIZE,NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NBOSIZE,
     +                NBOSIZE,NMOVE,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        B )
     +
     +
         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                COLMAP,
     +                IVEC,
     +
     +                        BDSIZE )
     +
     +
         DO I = 1,BONDSIZE

            CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +                 ( NMOVE,
     +                   NMOVE,
     +                   NMOVE,
     +                   NMOVE,
     +                   COLMAP,
     +                   IVEC,
     +
     +                           BDATOM (1,I) )
     +
     +
         END DO
C
C
C             ...form now the Core and remaining Rydberg NBO orbitals.
C                The array H containing the atomic NHOs will be used
C                for temporary space and the atomic NHO in it will be
C                destroyed.
C
C
         CALL  NLO__FORM_ATOMIC_NONVALENCE_NBO
     +
     +              ( NBAS,MXNBA,
     +                NCB,NCA,
     +                ATNCB,ATCOFF,
     +                P,
     +                H (1,NHB+1),
     +                XMAT,
     +
     +                        W (1),
     +                        C (1,1))
     +
     +
         CALL  NLO__FORM_ATOMIC_NONVALENCE_NBO
     +
     +              ( NBAS,MXNBA,
     +                NRB,NRA,
     +                ATNRB,ATROFF,
     +                P,
     +                H (1,NHB+NCB+1),
     +                XMAT,
     +
     +                        W (NHB+NCB+1),
     +                        C (1,NHB+NCB+1))
     +
     +
C
C
C             ...generate NBO info for Lone-pairs, Empty-pairs and
C                Rydberg NBOs. 
C
C
         OFFLP = NCB
         OFFEP = OFFLP + NLB + NBB
         OFFRY = OFFEP + NEB + NAB

         NLA = 0

         IF (NLB.GT.0) THEN
             NLBA = 1
             NLBSUM = 0
             ATOLD = BDATOM (OFFLP+1,1)

             DO I = 2,NLB
                ATNEW = BDATOM (OFFLP+I,1)
                IF (ATNEW.EQ.ATOLD) THEN
                    NLBA = NLBA + 1
                ELSE
                    NLA = NLA + 1
                    ATNLB (NLA) = NLBA
                    ATLIDX (NLA) = ATOLD
                    ATLOFF (NLA) = NLBSUM
                    NLBSUM = NLBSUM + NLBA
                    NLBA = 1
                    ATOLD = ATNEW
                END IF
             END DO

             NLA = NLA + 1
             ATNLB (NLA) = NLBA
             ATLIDX (NLA) = ATOLD
             ATLOFF (NLA) = NLBSUM
         END IF

         NEA = 0

         IF (NEB.GT.0) THEN
             NEBA = 1
             NEBSUM = 0
             ATOLD = BDATOM (OFFEP+1,1)

             DO I = 2,NEB
                ATNEW = BDATOM (OFFEP+I,1)
                IF (ATNEW.EQ.ATOLD) THEN
                    NEBA = NEBA + 1
                ELSE
                    NEA = NEA + 1
                    ATNEB (NEA) = NEBA
                    ATEIDX (NEA) = ATOLD
                    ATEOFF (NEA) = NEBSUM
                    NEBSUM = NEBSUM + NEBA
                    NEBA = 1
                    ATOLD = ATNEW
                END IF
             END DO

             NEA = NEA + 1
             ATNEB (NEA) = NEBA
             ATEIDX (NEA) = ATOLD
             ATEOFF (NEA) = NEBSUM
         END IF

         NYA = 0

         IF (NYB.GT.0) THEN
             NYBA = 1
             NYBSUM = 0
             ATOLD = BDATOM (OFFRY+1,1)

             DO I = 2,NYB
                ATNEW = BDATOM (OFFRY+I,1)
                IF (ATNEW.EQ.ATOLD) THEN
                    NYBA = NYBA + 1
                ELSE
                    NYA = NYA + 1
                    ATNYB (NYA) = NYBA
                    ATYIDX (NYA) = ATOLD
                    ATYOFF (NYA) = NYBSUM
                    NYBSUM = NYBSUM + NYBA
                    NYBA = 1
                    ATOLD = ATNEW
                END IF
             END DO

             NYA = NYA + 1
             ATNYB (NYA) = NYBA
             ATYIDX (NYA) = ATOLD
             ATYOFF (NYA) = NYBSUM
         END IF
C
C
C             ...calculate the NBO occupation interaction orders. 
C
C
         DO NBO = 1,NBAS

            CALL  NLO__CALC_INTERACTION_ORDER
     +
     +                 ( NBAS,NBAS,
     +                   NBAS,NBAS,
     +                   NBAS,
     +                   NBAS,NBAS,
     +                   NBAS,
     +                   NBO,
     +                   C,P,C (1,NBO),
     +                   XVEC (1),
     +                   XVEC (NBAS+1),
     +
     +                           Q (NBO) )
     +
     +
         END DO

         CALL    MAT__PRINT_A_FLOAT_18_NOZEROS
     +
     +                ( 6,
     +                  ' NBO Interaction order vector (not sym)',
     +                  NBAS,1,
     +                  NBAS,1,
     +                  Q )
     +
     +
C
C
C             ...evaluate the NBO angular momentum table.
C
C
         DO LTYPE = 0,MXSHELL

            DO NBO = 1,NBAS
               BOND = NBOBD (NBO)
               NNHO = BDNBAS (BOND)
               ANGVAL = ZERO

               DO I = 1,NNHO
                  NHO = BDBAS (I,BOND)
                  CSQR = B (I,NBO) ** 2
                  ANGVAL = ANGVAL + CSQR * ANGNHO (NHO,LTYPE)
               END DO

               ANGNBO (NBO,LTYPE) = ANGVAL

            END DO
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
