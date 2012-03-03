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
         SUBROUTINE  NLO__SYMMETRIZE_NBO_ORBITALS
     +
     +                    ( NBAS,NATOM,NBOND,
     +                      NBOSIZE,ROTSIZE,
     +                      MXNBA,MXSHELL,
     +                      NCB,NLB,NBB,NEB,NAB,NYB,NRB,
     +                      NCA,NLA,NEA,NYA,NRA,
     +                      ATNCB,ATCIDX,ATCOFF,
     +                      ATNLB,ATLIDX,ATLOFF,
     +                      ATNEB,ATEIDX,ATEOFF,
     +                      ATNYB,ATYIDX,ATYOFF,
     +                      ATNRB,ATRIDX,ATROFF,
     +                      BASBEG,BASEND,
     +                      DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                      SYMMAP,SYMNBO,SYMCEN,
     +                      PLATONIC,
     +                      RING,NRING,RINGSZ,
     +                      NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                      PLATO,
     +                      BDNCEN,BDCEN,
     +                      ANGNHO,
     +                      HYB,ROT,
     +                      P,
     +                      DOMAT,
     +                      SHALF,
     +                      LOCAL,
     +                      SLTYPE,
     +                      IVEC,
     +                      XVEC,XMAT,
     +
     +                              NBOBD,
     +                              BDNBAS,BDBAS,
     +                              ANGNBO,
     +                              B,W,Q,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__SYMMETRIZE_NBO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine symmetrizes the sets of natural bond
C                orbitals (NBOs).
C
C
C                  Input:
C
C                    NBAS         =  total # of AOs in AO basis
C                    NATOM        =  total # of atomic centers
C                    NBOND        =  total # of bonds found in the NHB
C                                    part.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    ROTSIZE      =  maximum # of NHOs that will be
C                                    allowed to form NBO hybrids for
C                                    rotationally invariant atomic NBO
C                                    construction.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    MXSHELL      =  largest l-shell value
C                    NxB          =  total # of Core, Lone-pair, Bond,
C                                    Empty-pair, Antibond and Rydberg
C                                    NBOs determined (x=C,L,B,E,A,Y,R)
C                    NxA          =  total # of Core, Lone-pair,
C                                    Empty-pair and Rydberg atoms
C                                    (x=C,L,E,Y,R)
C                    ATNCB (A)    =  # of Core NBOs on core atom A.
C                    ATCIDX (A)   =  atomic index for core atom A.
C                    ATCOFF (A)   =  index offset for Core NBOs for
C                                    core atom A. This index is equal
C                                    to the total number of Core NBOs
C                                    on all core atoms preceeding
C                                    core atom A.
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
C                    ATNYB (A)    =  # of hybrid Rydberg NBOs on hybrid
C                                    Rydberg atom A.
C                    ATYIDX (A)   =  atomic index for hybrid Rydberg
C                                    atom A.
C                    ATYOFF (A)   =  index offset for hybrid Rydberg
C                                    NBOs for hybrid Rydberg atom A.
C                                    This index is equal to the total
C                                    number of hybrid Rydberg NBOs on
C                                    all hybrid Rydberg atoms preceeding
C                                    hybrid Rydberg atom A.
C                    ATNRB (A)    =  # of Rydberg NBOs on Rydberg atom A
C                    ATRIDX (A)   =  atomic index for Rydberg atom A.
C                    ATROFF (A)   =  index offset for Rydberg NBOs for
C                                    Rydberg atom A. This index is
C                                    equal to the total number of
C                                    Rydberg NBOs on all Rydberg atoms
C                                    preceeding Rydberg atom A.
C                    BASBEG (A)   =  first basis index number for atom A
C                    BASEND (A)   =  last basis index number for atom A
C                    DSYMACC      =  symmetry accuracy for distance
C                                    order matrix values
C                    LSYMACC      =  symmetry accuracy for orbital
C                                    localization contents
C                    PSYMACC      =  symmetry accuracy for occupation
C                                    matrix values (including weights)
C                    QSYMACC      =  symmetry accuracy for orbital
C                                    interaction order values (sum of
C                                    squares of offdiagonal occupation
C                                    matrix elements for one orbital)
C                    SYMMAP (A,B) =  has been set equal to 1 if atoms
C                                    A and B were found to be symmetry
C                                    related at the current stage of
C                                    calculation. A value of 0 indicates
C                                    no symmetry relation.
C                    SYMNBO       =  integer vector which will be used
C                                    to indicate which NBO is already
C                                    considered to be symmetry adapted
C                                    at the different stages of
C                                    symmetrization.
C                    SYMCEN       =  array that will be used to store
C                                    symmetry related atomic indices.
C                    PLATONIC     =  is true, if the complete set of
C                                    atomic centers possess one of the
C                                    possible 5 platonic symmetries.
C                    RING (A,B)   =  integer array containing all atomic
C                                    indices A of those symmetry related
C                                    sets that have equal distance from
C                                    center B and are of maximum size
C                    NRING (B)    =  integer vector containing the # of
C                                    sets of maximum size related to
C                                    center B.
C                    RINGSZ (B)   =  integer vector containing the
C                                    maximum size of the above found
C                                    sets related to center B.
C                    NTETRA       =  contains the # of tetrahedrally
C                                    related atomic center sets (size 4)
C                    NCUBE        =  contains the # of cubically
C                                    related atomic center sets (size 8)
C                    NOCTA        =  contains the # of octahedrally
C                                    related atomic center sets (size 6)
C                    NDODECA      =  contains the # of dodecahedrally
C                                    related atomic center sets
C                                    (size 20)
C                    NICOSA       =  contains the # of icosahedrally
C                                    related atomic center sets
C                                    (size 12)
C                    PLATO (A,P)  =  integer array containing all atomic
C                                    indices A for specific platonic
C                                    symmetry types as characterized by
C                                    the index P: P = 1 (four atoms
C                                    tetrahedrally arranged), P = 2
C                                    (eight atoms arranged in a cube),
C                                    P = 3 (six octahedral atoms), P = 4
C                                    (twenty atoms in a dodecahedron),
C                                    P = 5 (twelve icosahedral atoms)
C                    BDNCEN (J)   =  # of atomic centers for J-th bond.
C                                    The bond order is NHB-bonds/NCB/
C                                    NRB. Note that the # of NHB-bonds
C                                    might be < NHB size.
C                    BDCEN (I,J)  =  I-th atomic center index for J-th
C                                    bond. The bond order is NHB-bonds/
C                                    NCB/NRB. Note that the # of
C                                    NHB-bonds might be < NHB size.
C                    ANGNHO (J,L) =  NBAS x (MXSHELL+1) matrix
C                                    containing the sum of the square
C                                    of the NAO coefficients per L-type
C                                    angular momentum for the J-th NHO.
C                                    The NHO order is NHB/NCB/NRB.
C                    HYB          =  will contain the NBO hybridization
C                                    matrix.
C                    ROT          =  will contain the NBO rotation
C                                    matrix.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    DOMAT (A,B)  =  distance order matrix containing
C                                    'distance' info between atoms A
C                                    and B. Will be used when attempting
C                                    NBO symmetrization.
C                    SHALF        =  full NBAS x NBAS square root of
C                                    the overlap matrix in AO basis.
C                    LOCAL        =  will hold NBO atomic localization
C                                    contents.
C                    SLTYPE (J)   =  will hold NBO smallest angular
C                                    momentum component.
C                    IVEC         =  int scratch array of vector type.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    NBOBD (I)    =  initial (before symmetrization)
C                                    bond index number for the I-th NBO.
C                                    This array is the handle for
C                                    accessing and modifying info of
C                                    the bonds sitting in the arrays
C                                    BDNCEN,BDCEN,BDBAS and BDOCC. The
C                                    NBOBD elements are in NCB/NLB/NBB/
C                                    NEB/NAB/NYB/NRB order.
C                    BDNBAS (J)   =  initial (before symmetrization)
C                                    # of basis functions (NHOs) for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  initial (before symmetrization)
C                                    I-th global basis (NHO) index for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    ANGNBO (J,L) =  initial (before symmetrization)
C                                    NBAS x (MXSHELL+1) matrix containing
C                                    the sum of the square of the NAO
C                                    coefficients per L-type angular
C                                    momentum for the J-th NBO. The NBO
C                                    order is NCB/NLB/NBB/NEB/NAB/NYB/
C                                    NRB.
C                    B (I,J)      =  initial (before symmetrization)
C                                    NBOSIZE x NBAS matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO. The NBO
C                                    column index is in NCB/NLB/NBB/NEB/
C                                    NAB/NYB/NRB order.
C                    W            =  initial (before symmetrization)
C                                    NBO weight vector with elements in
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB order.
C                    Q            =  initial (before symmetrization)
C                                    NBO occupation interaction order
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
C                    C            =  initial (before symmetrization)
C                                    full NBO coefficient matrix
C                                    in AO basis with columns in
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB order.
C
C
C                  Output:
C
C                    NBOBD (I)    =  final (after symmetrization)
C                                    bond index number for the I-th NBO.
C                                    This array is the handle for
C                                    accessing and modifying info of
C                                    the bonds sitting in the arrays
C                                    BDNCEN,BDCEN,BDBAS and BDOCC. The
C                                    NBOBD elements are in NCB/NLB/NBB/
C                                    NEB/NAB/NYB/NRB order.
C                    BDNBAS (J)   =  final (after symmetrization)
C                                    # of basis functions (NHOs) for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  final (after symmetrization)
C                                    I-th global basis (NHO) index for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    ANGNBO (J,L) =  final (after symmetrization)
C                                    NBAS x (MXSHELL+1) matrix containing
C                                    the sum of the square of the NAO
C                                    coefficients per L-type angular
C                                    momentum for the J-th NBO. The NBO
C                                    order is NCB/NLB/NBB/NEB/NAB/NYB/
C                                    NRB.
C                    B (I,J)      =  final (after symmetrization)
C                                    NBOSIZE x NBAS matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO. The NBO
C                                    column index is in NCB/NLB/NBB/NEB/
C                                    NAB/NYB/NRB order.
C                    W            =  final (after symmetrization)
C                                    NBO weight vector with elements in
C                                    NCB/NLB/NBB/NEB/NAB/NYB/NRB order.
C                    Q            =  final (after symmetrization)
C                                    NBO occupation interaction order
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
C                    C            =  final (after symmetrization)
C                                    full NBO coefficient matrix
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

         LOGICAL     INITSYM,ATOMSYM,AXISSYM,HIGHSYM
         LOGICAL     PLATONIC
         LOGICAL     SKIP
         LOGICAL     TOTAL

         INTEGER     BOND
         INTEGER     I
         INTEGER     LTYPE
         INTEGER     MXNBA,MXSHELL
         INTEGER     NBAS,NATOM
         INTEGER     NBO
         INTEGER     NBOND
         INTEGER     NBOSIZE,ROTSIZE
         INTEGER     NCEN
         INTEGER     NTETRA,NCUBE,NOCTA,NDODECA,NICOSA
         INTEGER     NCB,NCA,NRB,NRA
         INTEGER     NHO,NNHO
         INTEGER     NLB,NLA,NBB,NEB,NEA,NAB,NYB,NYA
         INTEGER     OFFCR,OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR
         INTEGER     SYMSTEP

         INTEGER     ATNCB   (1:NCA   )
         INTEGER     ATNRB   (1:NRA   )
         INTEGER     ATNLB   (1:NLA   )
         INTEGER     ATNEB   (1:NEA   )
         INTEGER     ATNYB   (1:NYA   )
         INTEGER     ATCIDX  (1:NCA   )
         INTEGER     ATRIDX  (1:NRA   )
         INTEGER     ATLIDX  (1:NLA   )
         INTEGER     ATEIDX  (1:NEA   )
         INTEGER     ATYIDX  (1:NYA   )
         INTEGER     ATCOFF  (1:NCA   )
         INTEGER     ATROFF  (1:NRA   )
         INTEGER     ATLOFF  (1:NLA   )
         INTEGER     ATEOFF  (1:NEA   )
         INTEGER     ATYOFF  (1:NYA   )
         INTEGER     BASBEG  (1:NATOM )
         INTEGER     BASEND  (1:NATOM )
         INTEGER     BDNBAS  (1:NBAS  )
         INTEGER     BDNCEN  (1:NBAS  )
         INTEGER     COLMAP  (1:NBAS  )
         INTEGER     IVEC    (1:3*NBAS)
         INTEGER     NBOBD   (1:NBAS  )
         INTEGER     NRING   (1:NATOM )
         INTEGER     RINGSZ  (1:NATOM )
         INTEGER     SLTYPE  (1:NBAS  )
         INTEGER     SYMCEN  (1:NATOM )
         INTEGER     SYMNBO  (1:NBAS  )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS )
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS )
         INTEGER     PLATO   (1:NATOM  ,1:5    )
         INTEGER     RING    (1:NATOM  ,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM  ,1:NATOM)

         DOUBLE PRECISION  ANGVAL
         DOUBLE PRECISION  CSQR
         DOUBLE PRECISION  DSYMACC,LSYMACC,PSYMACC,QSYMACC
         DOUBLE PRECISION  SYMCRIT
         DOUBLE PRECISION  ZERO,ONE,TINY

         DOUBLE PRECISION  LOCAL (1:NBAS     )
         DOUBLE PRECISION  Q     (1:NBAS     )
         DOUBLE PRECISION  W     (1:NBAS     )
         DOUBLE PRECISION  XVEC  (1:NBAS+NBAS)

         DOUBLE PRECISION  ANGNBO (1:NBAS   ,0:MXSHELL)
         DOUBLE PRECISION  ANGNHO (1:NBAS   ,0:MXSHELL)
         DOUBLE PRECISION  B      (1:NBOSIZE,1:NBAS   )
         DOUBLE PRECISION  C      (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  P      (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  HYB    (1:ROTSIZE,1:ROTSIZE)
         DOUBLE PRECISION  ROT    (1:ROTSIZE,1:ROTSIZE)
         DOUBLE PRECISION  SHALF  (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  DOMAT  (1:NATOM  ,1:NATOM  )
         DOUBLE PRECISION  XMAT   (1:NBAS   ,1:NBAS   )

         PARAMETER  (ZERO  = 0.D0 )
         PARAMETER  (ONE   = 1.D0 )
         PARAMETER  (TINY  = 1.D-6)
C
C
C------------------------------------------------------------------------
C
C
C             ...attempt symmetrization (if necessary). 
C
C
         SKIP  =        (NCB.EQ.0)
     +            .AND. (NLB.EQ.0)
     +            .AND. (NEB.EQ.0)
     +            .AND. (NYB.EQ.0)
     +            .AND. (NRB.EQ.0)

         IF (SKIP) THEN
             RETURN
         END IF

         OFFCR = 0
         OFFLP = NCB
         OFFBD = OFFLP + NLB
         OFFEP = OFFBD + NBB
         OFFAB = OFFEP + NEB
         OFFRY = OFFAB + NAB
         OFFRR = OFFRY + NYB
C
C
C             ...determine all NBO total atomic localization contents.
C                These will assist together with the NBO interaction
C                orders in identifying symmetry equivalent NBOs.
C
C
         TOTAL = .TRUE.

         DO NBO = 1,NBAS

            BOND = NBOBD (NBO)
            NCEN = BDNCEN (BOND)

            CALL  NLO__EXTRACT_ATOMIC_CONTENT
     +
     +                 ( NBAS,NATOM,
     +                   NCEN,
     +                   BDCEN (1,BOND),
     +                   BASBEG,BASEND,
     +                   TOTAL,
     +                   C (1,NBO),
     +                   SHALF,
     +
     +                            XVEC (1) )
     +
     +
            LOCAL (NBO) = XVEC (1)

         END DO

         CALL    MAT__PRINT_A_FLOAT_18_NOZEROS
     +
     +                ( 6,
     +                  ' Locality vector (for symmetry use) ',
     +                  NBAS,1,
     +                  NBAS,1,
     +                  LOCAL )
     +
     +
C
C
C             ...for symmetrization we need to temporarily evaluate
C                the NBO smallest angular momentum vector, which will
C                possibly be changed during the symmetrization process.
C                Outermost loop is over the angular momentum type
C                in decreasing order!
C
C
         DO LTYPE = MXSHELL,0,-1
            DO NBO = 1,NBAS
               ANGVAL = ANGNBO (NBO,LTYPE)
               IF (ANGVAL.GT.TINY) THEN
                   SLTYPE (NBO) = LTYPE
               END IF
            END DO
         END DO
C
C
C             ...set initital values for array SYMNBO. Note, that
C                we use only bonds and anti-bonds of size 2 for
C                symmetrization. Larger sized bonds and anti-bonds
C                are not yet symmetry adapted and cannot be used
C                at this initial stage.
C
C
         DO I = 1,NBAS
            SYMNBO (I) = 0
         END DO

         DO I = 1,NBB
            NBO = OFFBD + I
            BOND = NBOBD (NBO)
            NCEN = BDNCEN (BOND)
            IF (NCEN.EQ.2) THEN
                SYMNBO (NBO) = 2
            END IF
         END DO

         DO I = 1,NAB
            NBO = OFFAB + I
            BOND = NBOBD (NBO)
            NCEN = BDNCEN (BOND)
            IF (NCEN.EQ.2) THEN
                SYMNBO (NBO) = 2
            END IF
         END DO
C
C
C             ...loop over all the symmetrization steps. The first
C                pass (steps 1-4) is done with a relatively high
C                symmetry interaction criterion, to avoid rotations
C                of degenerate NBOs using only possible computational
C                noise background. For the second pass this criterion
C                is set equal to the accuracy of the occupation matrix
C                to achieve symmetrization with all what we have at
C                hand by then.
C
C
         DO SYMSTEP = 1,6

            INITSYM = SYMSTEP.EQ.1
            ATOMSYM = SYMSTEP.EQ.2 .OR. SYMSTEP.EQ.5
            AXISSYM = SYMSTEP.EQ.3 .OR. SYMSTEP.EQ.6
            HIGHSYM = SYMSTEP.EQ.4

            IF (SYMSTEP.LE.4) THEN
                SYMCRIT = DSQRT (PSYMACC)
            ELSE
                SYMCRIT = PSYMACC
            END IF
C
C
C             ...Core NBO symmetrization (if necessary). 
C
C
            IF (NCB.GT.0) THEN

                CALL  NLO__SYMMETRIZE_ATOMIC_NBO
     +
     +                     ( NBAS,NATOM,
     +                       NBOSIZE,ROTSIZE,
     +                       MXNBA,
     +                       NCB,NCA,
     +                       OFFCR,
     +                       OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR,
     +                       ATNCB,ATCIDX,ATCOFF,IVEC (1),
     +                       BASBEG,BASEND,
     +                       COLMAP,SYMMAP,
     +                       IVEC (NBAS+1),
     +                       INITSYM,ATOMSYM,AXISSYM,HIGHSYM,
     +                       DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                       PLATONIC,
     +                       RING,NRING,RINGSZ,
     +                       NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                       PLATO,
     +                       SYMCEN,
     +                       BDNCEN,BDCEN,
     +                       HYB,
     +                       ROT,IVEC (2*NBAS+1),
     +                       SYMCRIT,
     +                       P,
     +                       DOMAT,
     +                       SHALF,
     +                       XVEC,
     +                       XMAT,
     +
     +                               SYMNBO,
     +                               SLTYPE,
     +                               NBOBD,
     +                               BDNBAS,
     +                               BDBAS,
     +                               LOCAL,
     +                               B,W,Q,C )
     +
     +
            END IF
C
C
C             ...Lone-pair NBO symmetrization (if necessary). 
C
C
            IF (NLB.GT.0) THEN

                CALL  NLO__SYMMETRIZE_ATOMIC_NBO
     +
     +                     ( NBAS,NATOM,
     +                       NBOSIZE,ROTSIZE,
     +                       MXNBA,
     +                       NLB,NLA,
     +                       OFFLP,
     +                       OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR,
     +                       ATNLB,ATLIDX,ATLOFF,IVEC (1),
     +                       BASBEG,BASEND,
     +                       COLMAP,SYMMAP,
     +                       IVEC (NBAS+1),
     +                       INITSYM,ATOMSYM,AXISSYM,HIGHSYM,
     +                       DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                       PLATONIC,
     +                       RING,NRING,RINGSZ,
     +                       NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                       PLATO,
     +                       SYMCEN,
     +                       BDNCEN,BDCEN,
     +                       HYB,
     +                       ROT,IVEC (2*NBAS+1),
     +                       SYMCRIT,
     +                       P,
     +                       DOMAT,
     +                       SHALF,
     +                       XVEC,
     +                       XMAT,
     +
     +                               SYMNBO,
     +                               SLTYPE,
     +                               NBOBD,
     +                               BDNBAS,
     +                               BDBAS,
     +                               LOCAL,
     +                               B,W,Q,C )
     +
     +
            END IF
C
C
C             ...Empty-pair NBO symmetrization (if necessary). 
C
C
            IF (NEB.GT.0) THEN

                CALL  NLO__SYMMETRIZE_ATOMIC_NBO
     +
     +                     ( NBAS,NATOM,
     +                       NBOSIZE,ROTSIZE,
     +                       MXNBA,
     +                       NEB,NEA,
     +                       OFFEP,
     +                       OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR,
     +                       ATNEB,ATEIDX,ATEOFF,IVEC (1),
     +                       BASBEG,BASEND,
     +                       COLMAP,SYMMAP,
     +                       IVEC (NBAS+1),
     +                       INITSYM,ATOMSYM,AXISSYM,HIGHSYM,
     +                       DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                       PLATONIC,
     +                       RING,NRING,RINGSZ,
     +                       NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                       PLATO,
     +                       SYMCEN,
     +                       BDNCEN,BDCEN,
     +                       HYB,
     +                       ROT,IVEC (2*NBAS+1),
     +                       SYMCRIT,
     +                       P,
     +                       DOMAT,
     +                       SHALF,
     +                       XVEC,
     +                       XMAT,
     +
     +                               SYMNBO,
     +                               SLTYPE,
     +                               NBOBD,
     +                               BDNBAS,
     +                               BDBAS,
     +                               LOCAL,
     +                               B,W,Q,C )
     +
     +
            END IF
C
C
C             ...Rydberg bond NBO symmetrization (if necessary). 
C
C
            IF (NYB.GT.0) THEN

                CALL  NLO__SYMMETRIZE_ATOMIC_NBO
     +
     +                     ( NBAS,NATOM,
     +                       NBOSIZE,ROTSIZE,
     +                       MXNBA,
     +                       NYB,NYA,
     +                       OFFRY,
     +                       OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR,
     +                       ATNYB,ATYIDX,ATYOFF,IVEC (1),
     +                       BASBEG,BASEND,
     +                       COLMAP,SYMMAP,
     +                       IVEC (NBAS+1),
     +                       INITSYM,ATOMSYM,AXISSYM,HIGHSYM,
     +                       DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                       PLATONIC,
     +                       RING,NRING,RINGSZ,
     +                       NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                       PLATO,
     +                       SYMCEN,
     +                       BDNCEN,BDCEN,
     +                       HYB,
     +                       ROT,IVEC (2*NBAS+1),
     +                       SYMCRIT,
     +                       P,
     +                       DOMAT,
     +                       SHALF,
     +                       XVEC,
     +                       XMAT,
     +
     +                               SYMNBO,
     +                               SLTYPE,
     +                               NBOBD,
     +                               BDNBAS,
     +                               BDBAS,
     +                               LOCAL,
     +                               B,W,Q,C )
     +
     +
            END IF
C
C
C             ...Rydberg NBO symmetrization (if necessary). 
C
C
            IF (NRB.GT.0) THEN

                CALL  NLO__SYMMETRIZE_ATOMIC_NBO
     +
     +                     ( NBAS,NATOM,
     +                       NBOSIZE,ROTSIZE,
     +                       MXNBA,
     +                       NRB,NRA,
     +                       OFFRR,
     +                       OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR,
     +                       ATNRB,ATRIDX,ATROFF,IVEC (1),
     +                       BASBEG,BASEND,
     +                       COLMAP,SYMMAP,
     +                       IVEC (NBAS+1),
     +                       INITSYM,ATOMSYM,AXISSYM,HIGHSYM,
     +                       DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                       PLATONIC,
     +                       RING,NRING,RINGSZ,
     +                       NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                       PLATO,
     +                       SYMCEN,
     +                       BDNCEN,BDCEN,
     +                       HYB,
     +                       ROT,IVEC (2*NBAS+1),
     +                       SYMCRIT,
     +                       P,
     +                       DOMAT,
     +                       SHALF,
     +                       XVEC,
     +                       XMAT,
     +
     +                               SYMNBO,
     +                               SLTYPE,
     +                               NBOBD,
     +                               BDNBAS,
     +                               BDBAS,
     +                               LOCAL,
     +                               B,W,Q,C )
     +
     +
            END IF

         END DO
C
C
C             ...reevaluate the NBO angular momentum table.
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
C             ...printout for checking. 
C
C
         CALL    MAT__PRINT_V_INTEGER_NOZEROS
     +
     +            ( 1,
     +              ' BDNBAS vector (after symmetrize) ',
     +              NBAS,
     +              NBAS,
     +              BDNBAS )
     +
     +
         CALL    MAT__PRINT_V_INTEGER_NOZEROS
     +
     +            ( 1,
     +              ' BDNCEN vector (after symmetrize) ',
     +              NBAS,
     +              NBAS,
     +              BDNCEN )
     +
     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +            ( 1,
     +              ' BDCEN array (after symmetrize) ',
     +              NBOSIZE,NBAS,
     +              NBOSIZE,NBAS,
     +              BDCEN )
     +
     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +            ( 1,
     +              ' BDBAS array (after symmetrize) ',
     +              NBOSIZE,NBAS,
     +              NBOSIZE,NBAS,
     +              BDBAS )
     +
     +
         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
     +
     +            ( 1,
     +              ' B matrix (after symmetrize) ',
     +              NBOSIZE,NBAS,
     +              NBOSIZE,NBAS,
     +              B )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
