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
         SUBROUTINE  NLO__ROTATE_ATOM_TGROUP_NBO
     +
     +                    ( NBAS,NATOM,NDEG,
     +                      NBOSIZE,ROTSIZE,
     +                      NXBA,
     +                      OFFX,OFFA,
     +                      OFFLP,OFFBD,OFFEP,
     +                      OFFAB,OFFRY,OFFRR,
     +                      CENTER,
     +                      NTETRA,
     +                      PLATO,
     +                      SYMCEN,
     +                      LSYMACC,QSYMACC,
     +                      FIRST,
     +                      IDXDEG,
     +                      NBOBD,
     +                      BDNCEN,BDCEN,
     +                      HYB,
     +                      ROT,ROTIDX,
     +                      SYMCRIT,
     +                      P,
     +                      LOCAL,
     +                      W,Q,
     +                      XVEC,
     +                      XMAT,
     +
     +                              SYMNBO,
     +                              SLTYPE,
     +                              BDNBAS,
     +                              BDBAS,
     +                              B,
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ROTATE_ATOM_TGROUP_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine rotates a degenerate set of NBOs sitting
C                at the center of tetrahedral symmetry. The background
C                atomic sets used for performing this rotation are
C                chosen from the tetrahedrally related centers with
C                C3v site symmetry. The following Td -> C3v irrep
C                decomposition table shows what will be expected:
C
C
C                                Td  |    C3v
C                            ------------------
C                                A1  |    A1
C                                A2  |    A2
C                                 E  |     E
C                                T1  |  A2 + E
C                                T2  |  A1 + E
C
C
C                From this table we can construct a table, showing
C                which tetrahedral degenerate NBO sets can have (+)
C                or cannot have (-) interactions with the tetrahedral
C                (C3v site symmetry) background:
C
C
C
C                      C3v background | Td-irrep: E  T1  T2
C                    ----------------------------------------
C                      A1 (s,pz)      |           -   -   +
C                      A2 (f,g,h,...) |           -   +   -
C                    A1+E (sp2-hybrid)|           +   +   +
C
C
C
C                The table also shows what type of C3v irrep will occur
C                for specific atomic functions used on those sites.
C                For example, A2 irreps of C3v will only occur when
C                including f- or higher atomic functions.
C
C                Another table is also useful, showing decomposition
C                of atomic angular momentum levels into Td irreps:
C
C
C                        atomic level  |          Td irreps
C                      -----------------------------------------------
C                             s        |              A1
C                             p        |              T2
C                             d        |            E + T2
C                             f        |         A2 + T1 + T2
C                             g        |       A1 + E + T1 + T2
C                             h        |          E + T1 + 2T2
C                             i        |    A1 + A2 + E + T1 + 2T2
C
C
C                We have next to distinguish two cases:
C
C                  i) 2-fold degenerate NBOs
C
C                     These can be only of E type and can occur if the
C                     basis set on the center atom has at least d-type
C                     functions.
C
C                     Procedure for the E type NBOs:
C
C                       a) Rotate the E pair to one of the background
C                          sites, such that the first E NBO interaction
C                          with this site is maximized (C3v A1+E sp2
C                          hybrid site symmetry).
C
C                       b) Once appropriately rotated, form sd2
C                          hybrids to achieve octahedral symmetry.
C                           
C
C                 ii) 3-fold degenerate NBOs
C
C                     These can be of T1 and T2 type. For T2 we can
C                     use the tetrahedral A1 background, for T1 we
C                     can use either the A2 (rare!) or the A1+E
C                     sp2-hybrid tetrahedral background. The 1-norm
C                     of the interaction occupation matrix between
C                     the background and the 3-fold degenerate NBOs
C                     is maximized.
C                     
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    NDEG         =  # of degenerate NBOs. Can only
C                                    be equal to 2,3,4 or 5.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    ROTSIZE      =  maximum # of NHOs that will be
C                                    allowed to form NBO hybrids for
C                                    rotationally invariant atomic NBO
C                                    construction.
C                    NXBA         =  # of atomic x-type NBOs on the
C                                    central atom treated here.
C                    OFFX         =  absolute offset index for present
C                                    x-type NBOs. This index indicates
C                                    the totallity of NBOs before the
C                                    x-type NBO set.
C                    OFFA         =  local atomic offset index for
C                                    present atomic x-type NBOs. This
C                                    index is equal to the total number
C                                    of atomic x-type NBOs on all
C                                    x-type NBO atoms preceeding x-type
C                                    NBO atom A.
C                    OFFxy        =  absolute offset indices for x-type
C                                    NBOs. These indices indicate the
C                                    totallity of x-type NBOs before the
C                                    x-type NBO set (xy:LP=Lone-pair,
C                                    BD=Bond, EP=Empty-pair, AB=Anti-
C                                    bond, RY,RR=Rydberg types)
C                    CENTER       =  atomic index of the center on which
C                                    NBO rotation will be performed
C                    NTETRA       =  contains the # of tetrahedrally
C                                    related center sets (size 4)
C                    PLATO (A)    =  contains all atomic indices A for
C                                    all tetrahedrally arranged centers
C                                    in groups of 4.
C                    SYMCEN       =  array that will be used to store
C                                    symmetry related atomic indices.
C                    LSYMACC      =  symmetry accuracy for orbital
C                                    localization contents
C                    QSYMACC      =  symmetry accuracy for orbital
C                                    interaction order values (sum of
C                                    squares of offdiagonal occupation
C                                    matrix elements for one orbital)
C                    FIRST        =  index of 1st degenerate NBO inside
C                                    the atomic NXBA-dimensional NBO set.
C                    IDXDEG       =  array that will be used to store
C                                    the NBO indices that will be
C                                    rotated (for example, if sd2
C                                    hybrids will be constructed, this
C                                    array will indicate all the NBO
C                                    indices for the s and the two d
C                                    NBOs).
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    for the I-th atomic x-type NBO.
C                                    This array is the handle for
C                                    accessing and modifying info of
C                                    the NHO bonds sitting in the
C                                    arrays BDNCEN,BDCEN,BDBAS and
C                                    BDOCC.
C                    BDNCEN (J)   =  # of atomic centers for J-th bond.
C                                    The bond order is NHB-bonds/NCB/
C                                    NRB. Note that the # of NHB-bonds
C                                    might be < NHB size.
C                    BDCEN (I,J)  =  I-th atomic center index for J-th
C                                    bond. The bond order is NHB-bonds/
C                                    NCB/NRB. Note that the # of
C                                    NHB-bonds might be < NHB size.
C                    HYB          =  will contain the NBO hybridization
C                                    matrix.
C                    ROT          =  will contain the NBO rotation
C                                    matrix.
C                    ROTIDX       =  will contain the NBO rotation
C                                    indices.
C                    SYMCRIT      =  symmetry interaction criterion.
C                                    Any absolute occupation matrix
C                                    interaction element between
C                                    symmetrized NBOs and a degenerate
C                                    set of NBOs is assumed to be zero
C                                    if it falls below this value.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    LOCAL        =  atomic localization content for
C                                    all NBOs.
C                    W            =  complete NBO weight vector at
C                                    present stage of symmetrization.
C                    Q            =  complete NBO interaction order
C                                    vector at present stage of
C                                    symmetrization.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    SYMNBO       =  initial integer vector indicating
C                                    the different symmetry stages of
C                                    each NBO (=0 means no symmetry
C                                    adaptation, =1 means hybridization
C                                    was performed, =2 means symmetry
C                                    adapted in final form). This array
C                                    is extremely useful for deciding
C                                    which NBOs can be used for the
C                                    symmetry adaptation at different
C                                    stages. Note that only NBOs with
C                                    SYMNBO = 2 can be considered for
C                                    that purpose.
C                    SLTYPE       =  initial integer vector indicating
C                                    the smallest angular momentum
C                                    component present in the NBOs.
C                                    Useful for deciding which atomic
C                                    NBOs to use for hybrid formation.
C                                    Once an atomic NBO has been used
C                                    for hybrid formation, a -1 will
C                                    be placed at the corresponding
C                                    vector position and the NBO cannot
C                                    be used again for hybrid formation.
C                    BDNBAS (J)   =  initial # of basis functions (NHOs)
C                                    for J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  initial I-th global basis (NHO)
C                                    index for J-th bond. The bond order
C                                    is NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    B (I,J)      =  initial NBOSIZE x NBAS matrix
C                                    containing the newly J-th symmetry
C                                    adapted x-type NBO expansion
C                                    coefficients in terms of the I-th
C                                    atomic NHOs forming the J-th
C                                    symmetry adapted x-type NBO.
C                    C            =  complete NBAS x NBAS NBO
C                                    coefficient matrix in AO basis
C                                    before symmetrization.
C
C
C                  Output:
C
C                    SYMNBO       =  updated integer vector indicating
C                                    the different symmetry stages of
C                                    each NBO (=0 means no symmetry
C                                    adaptation, =1 means hybridization
C                                    was performed, =2 means symmetry
C                                    adapted in final form). This array
C                                    is extremely useful for deciding
C                                    which NBOs can be used for the
C                                    symmetry adaptation at different
C                                    stages. Note that only NBOs with
C                                    SYMNBO = 2 can be considered for
C                                    that purpose.
C                    SLTYPE       =  updated integer vector indicating
C                                    the smallest angular momentum
C                                    component present in the NBOs.
C                                    Useful for deciding which atomic
C                                    NBOs to use for hybrid formation.
C                                    Once an atomic NBO has been used
C                                    for hybrid formation, a -1 will
C                                    be placed at the corresponding
C                                    vector position and the NBO cannot
C                                    be used again for hybrid formation.
C                    BDNBAS (J)   =  modified # of basis functions
C                                    (NHOs) for J-th bond. The bond
C                                    order is NHB-bonds/NCB/NRB. Note
C                                    that the # of NHB-bonds might
C                                    be < NHB size.
C                    BDBAS (I,J)  =  modified I-th global basis (NHO)
C                                    index for J-th bond. The bond order
C                                    is NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    B (I,J)      =  modified NBOSIZE x NBAS matrix
C                                    containing the newly J-th symmetry
C                                    adapted x-type NBO expansion
C                                    coefficients in terms of the I-th
C                                    atomic NHOs forming the J-th
C                                    symmetry adapted x-type NBO.
C                    C            =  complete NBAS x NBAS NBO
C                                    coefficient matrix in AO basis
C                                    after symmetrization.
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

         LOGICAL     HYBRID
         LOGICAL     MAXIMZE
         LOGICAL     NORMLZE
         LOGICAL     ROT2DEG,ROT3DEG
         LOGICAL     SAVEC,SAVEP
         LOGICAL     ZEROSYM

         INTEGER     BOND1,BOND2,BOND3
         INTEGER     CENTER
         INTEGER     D,N,S
         INTEGER     D1,D2
         INTEGER     FIRST
         INTEGER     H1,H2,H3
         INTEGER     IDX
         INTEGER     LOBE
         INTEGER     NROT
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NBOSIZE,ROTSIZE
         INTEGER     NHO1,NHO2,NHO3
         INTEGER     NDEG
         INTEGER     NITER
         INTEGER     NLO__ANGULAR_MOMENTUM_NBO
         INTEGER     NSYM
         INTEGER     NTETRA
         INTEGER     NXBA
         INTEGER     OFFA,OFFX
         INTEGER     OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR
         INTEGER     SYMSIZE,SYMSETS
         INTEGER     T1,T2,T3
         INTEGER     TOTSYM

         INTEGER     BDNBAS  (1:NBAS )
         INTEGER     BDNCEN  (1:NBAS )
         INTEGER     IDXDEG  (1:NXBA )
         INTEGER     NBOBD   (1:NBAS )
         INTEGER     PLATO   (1:NATOM)
         INTEGER     ROTIDX  (1:NBAS )
         INTEGER     SLTYPE  (1:NBAS )
         INTEGER     SYMCEN  (1:NATOM)
         INTEGER     SYMNBO  (1:NBAS )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS )
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS )

         DOUBLE PRECISION  LSYMACC,QSYMACC
         DOUBLE PRECISION  R11,R12,R21,R22
         DOUBLE PRECISION  SYMCRIT
         DOUBLE PRECISION  X,Y,Z
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  LOCAL (1:NBAS)
         DOUBLE PRECISION  Q     (1:NBAS)
         DOUBLE PRECISION  W     (1:NBAS)
         DOUBLE PRECISION  XVEC  (1:NBAS+NBAS)

         DOUBLE PRECISION  B     (1:NBOSIZE,1:NBAS   )
         DOUBLE PRECISION  C     (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  HYB   (1:ROTSIZE,1:ROTSIZE)
         DOUBLE PRECISION  P     (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  ROT   (1:ROTSIZE,1:ROTSIZE)
         DOUBLE PRECISION  XMAT  (1:NBAS   ,1:ROTSIZE)

         PARAMETER  (ZERO = 0.D0)
         PARAMETER  (ONE  = 1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check symmetry consistency. Since this routine deals
C                with the tetrahedral NBOs, we know that deg = 2 or 3
C                always. Otherwise stop with message.
C
C
         IF (NDEG.NE.2 .AND. NDEG.NE.3) THEN
             WRITE (*,*) ' Degeneracy conflict for atomic tgroup case! '
             WRITE (*,*) ' NDEG = ',NDEG
             WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
             WRITE (1,*) ' Degeneracy conflict for atomic tgroup case! '
             WRITE (1,*) ' NDEG = ',NDEG
             WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
             STOP
         END IF
C
C
C             ...check dimensions.
C
C
         IF (NDEG.GT.ROTSIZE) THEN
             WRITE (*,*) ' Degeneracy size too large! '
             WRITE (*,*) ' NDEG,ROTSIZE = ',NDEG,ROTSIZE
             WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
             WRITE (1,*) ' Degeneracy size too large! '
             WRITE (1,*) ' NDEG,ROTSIZE = ',NDEG,ROTSIZE
             WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
             STOP
         END IF
C
C
C             ...set some data.
C
C
         SAVEC = .TRUE.
         SAVEP = .TRUE.
C
C
C             ...handle 2-fold degeneracy, if any.
C
C
         IF (NDEG.EQ.2) THEN

             D = OFFX + OFFA + FIRST

             CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                  ( NBAS,NBAS,
     +                    NBAS,NBAS,
     +                    NBAS,2,
     +                    NBAS,2,
     +                    NBAS,
     +                    NBAS,NBAS,2,
     +                    0,0,
     +                    SAVEC,SAVEP,SAVEC,
     +                    C,P,C (1,D),
     +                    XVEC,
     +
     +                             XMAT )
     +
     +
             SYMSIZE = 4
             SYMSETS = NTETRA

             NSYM = SYMSIZE * SYMSETS

             DO N = 1,NSYM
                SYMCEN (N) = PLATO (N)
             END DO

             CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                  ( NBAS,NATOM,2,
     +                    NBOSIZE,
     +                    OFFLP,OFFBD,OFFEP,
     +                    OFFAB,OFFRY,OFFRR,
     +                    CENTER,
     +                    SYMSIZE,SYMSETS,SYMCEN,
     +                    SYMNBO,
     +                    QSYMACC,LSYMACC,
     +                    BDNCEN,BDCEN,
     +                    NBOBD,
     +                    SYMCRIT,
     +                    XMAT,
     +                    LOCAL,
     +                    Q,
     +
     +                             ZEROSYM,
     +                             NROT,
     +                             ROTIDX )
     +
     +
             ROT2DEG = .NOT.ZEROSYM

             IF (ROT2DEG) THEN
                 IF (NROT.EQ.0) THEN
                     WRITE (*,*) ' No 2-deg tgroup NBO rot set! '
                     WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
                     WRITE (1,*) ' No 2-deg tgroup NBO rot set! '
                     WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
                     STOP
                 END IF
             ELSE
                 WRITE (*,*) ' Could not rotate 2-deg tgroup NBOs! '
                 WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
                 WRITE (1,*) ' Could not rotate 2-deg tgroup NBOs! '
                 WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
                 RETURN
             END IF
C
C
C             ...up to now, the 2-fold basis forming the E-irreps
C                is mixed. For example, the 2z*z-x*x-y*y and x*x-y*y
C                d-functions are randomly mixed, thus each of the
C                d-functions possesses variable amounts of the z*z
C                component. In order to form clean sd2 hybrids, we
C                have to unmix the d-functions such that the first
C                one is identified as the one depleted of any z*z
C                component. This is achieved by rotating both
C                d-functions such that the first has no interaction
C                with the first background point. This first d-function
C                is then taken as the pure x*x-y*y type function
C                and passed to the sd2 hybrid formation routine.
C
C
             LOBE = ROTIDX (1)
             X    = XMAT (LOBE,1)
             Y    = XMAT (LOBE,2)

             MAXIMZE = .FALSE.

             CALL  NLO__FIND_2DEG_AXIAL_ROTMAT
     +
     +                  ( X,Y,
     +                    MAXIMZE,
     +
     +                         R11,R12,
     +                         R21,R22 )
     +
     +
             DO N = 1,NBAS
                X = R11 * C (N,D) + R21 * C (N,D+1)
                Y = R12 * C (N,D) + R22 * C (N,D+1)
                C (N,D  ) = X
                C (N,D+1) = Y
             END DO
C
C
C             ...if possible, form the symmetric octahedral sd2
C                hybrids.
C
C
             HYBRID = .FALSE.

             TOTSYM = NLO__ANGULAR_MOMENTUM_NBO
     +
     +                     ( NXBA,
     +                       SLTYPE (OFFX+OFFA+1),
     +                       FIRST,FIRST+1,
     +                       0,
     +                       1,
     +                       QSYMACC,LSYMACC,
     +                       W (OFFX+OFFA+1),
     +                       Q (OFFX+OFFA+1),
     +                       LOCAL (OFFX+OFFA+1) )
     +
     +
             HYBRID = TOTSYM .NE. 0

             IF (HYBRID) THEN

                 IF (ROTSIZE.LT.3) THEN
                     WRITE (*,*) ' Cannot form sd2 hybrids! '
                     WRITE (*,*) ' ROTSIZE = ',ROTSIZE
                     WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
                     WRITE (1,*) ' Cannot form sd2 hybrids! '
                     WRITE (1,*) ' ROTSIZE = ',ROTSIZE
                     WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
                     STOP
                 END IF

                 S = OFFX + OFFA + TOTSYM

                 ROT (1,1) = ZERO
                 ROT (2,1) = ZERO
                 ROT (3,1) = ZERO
                 ROT (1,2) = ZERO
                 ROT (2,2) = R11
                 ROT (3,2) = R21
                 ROT (1,3) = ZERO
                 ROT (2,3) = R12
                 ROT (3,3) = R22

                 CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                      ( NBAS,NBAS,
     +                        NBAS,NBAS,
     +                        NBAS,1,
     +                        NBAS,1,
     +                        NBAS,
     +                        NBAS,NBAS,1,
     +                        0,0,
     +                        SAVEC,SAVEP,SAVEC,
     +                        C,P,C (1,S),
     +                        XVEC,
     +
     +                                 XMAT )
     +
     +
                 X = XMAT (LOBE,1)

                 IF (X.LT.ZERO) THEN
                     DO N = 1,NBAS
                        C (N,S) =  - C (N,S)
                     END DO
                     ROT (1,1) = - ONE
                 ELSE
                     ROT (1,1) = ONE
                 END IF

                 CALL  NLO__FORM_SD2_NBO
     +
     +                      ( NBAS,
     +
     +                               HYB (1,1),HYB (1,2),HYB (1,3),
     +                               HYB (2,1),HYB (2,2),HYB (2,3),
     +                               HYB (3,1),HYB (3,2),HYB (3,3),
     +                               C (1,S),
     +                               C (1,D) )
     +
     +
                 H1 = S
                 H2 = D
                 H3 = D + 1

                 CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                      ( NBAS,NBAS,
     +                        NBAS,NBAS,
     +                        NBAS,1,
     +                        NBAS,1,
     +                        NBAS,
     +                        NBAS,NBAS,1,
     +                        0,0,
     +                        SAVEC,SAVEP,SAVEC,
     +                        C,P,C (1,S),
     +                        XVEC,
     +
     +                                 XMAT (1,1) )
     +
     +
                 CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                      ( NBAS,NBAS,
     +                        NBAS,NBAS,
     +                        NBAS,2,
     +                        NBAS,2,
     +                        NBAS,
     +                        NBAS,NBAS,2,
     +                        0,0,
     +                        SAVEC,SAVEP,SAVEC,
     +                        C,P,C (1,D),
     +                        XVEC,
     +
     +                                 XMAT (1,2) )
     +
     +
                 DO N = 1,NROT
                    IDX = ROTIDX (N)
                    XMAT (N,1) = XMAT (IDX,1)
                    XMAT (N,2) = XMAT (IDX,2)
                    XMAT (N,3) = XMAT (IDX,3)
                 END DO

                 CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                      ( 6,
     +                        ' P for sd2 orientation (atom,tgroup) ',
     +                        NBAS,3,
     +                        NROT,3,
     +                        XMAT )
     +
     +
C
C
C             ...form the hybrid NBO coefficient matrix in NHO basis and
C                the matrices containing details about the new NBOs.
C                Mark the new hybrid NBOs as symmetry adapted.
C
C
                 B (1,H1) =   ROT (1,1) * HYB (1,1)
     +                      + ROT (1,2) * HYB (2,1)
     +                      + ROT (1,3) * HYB (3,1)
                 B (2,H1) =   ROT (2,1) * HYB (1,1)
     +                      + ROT (2,2) * HYB (2,1)
     +                      + ROT (2,3) * HYB (3,1)
                 B (3,H1) =   ROT (3,1) * HYB (1,1)
     +                      + ROT (3,2) * HYB (2,1)
     +                      + ROT (3,3) * HYB (3,1)
                 B (1,H2) =   ROT (1,1) * HYB (1,2)
     +                      + ROT (1,2) * HYB (2,2)
     +                      + ROT (1,3) * HYB (3,2)
                 B (2,H2) =   ROT (2,1) * HYB (1,2)
     +                      + ROT (2,2) * HYB (2,2)
     +                      + ROT (2,3) * HYB (3,2)
                 B (3,H2) =   ROT (3,1) * HYB (1,2)
     +                      + ROT (3,2) * HYB (2,2)
     +                      + ROT (3,3) * HYB (3,2)
                 B (1,H3) =   ROT (1,1) * HYB (1,3)
     +                      + ROT (1,2) * HYB (2,3)
     +                      + ROT (1,3) * HYB (3,3)
                 B (2,H3) =   ROT (2,1) * HYB (1,3)
     +                      + ROT (2,2) * HYB (2,3)
     +                      + ROT (2,3) * HYB (3,3)
                 B (3,H3) =   ROT (3,1) * HYB (1,3)
     +                      + ROT (3,2) * HYB (2,3)
     +                      + ROT (3,3) * HYB (3,3)

                 BOND1 = NBOBD (H1)
                 BOND2 = NBOBD (H2)
                 BOND3 = NBOBD (H3)

                 NHO1 = BDBAS (1,BOND1)
                 NHO2 = BDBAS (1,BOND2)
                 NHO3 = BDBAS (1,BOND3)

                 BDNBAS (BOND1) = 3
                 BDNBAS (BOND2) = 3
                 BDNBAS (BOND3) = 3

                 BDBAS (1,BOND1) = NHO1
                 BDBAS (2,BOND1) = NHO2
                 BDBAS (3,BOND1) = NHO3
                 BDBAS (1,BOND2) = NHO1
                 BDBAS (2,BOND2) = NHO2
                 BDBAS (3,BOND2) = NHO3
                 BDBAS (1,BOND3) = NHO1
                 BDBAS (2,BOND3) = NHO2
                 BDBAS (3,BOND3) = NHO3

                 SYMNBO (H1) = 2
                 SYMNBO (H2) = 2
                 SYMNBO (H3) = 2

                 SLTYPE (H1) = -1
                 SLTYPE (H2) = -1
                 SLTYPE (H3) = -1

             ELSE

                 WRITE (*,*) ' No s function for sd2 hybrid! '
                 WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
                 WRITE (*,*) ' T-group symmetry reduction will occur! '
                 WRITE (1,*) ' No s function for sd2 hybrid! '
                 WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
                 WRITE (1,*) ' T-group symmetry reduction will occur! '

                 D1 = D
                 D2 = D + 1
C
C
C             ...form the 2-deg NBO coefficient matrix in NHO basis and
C                the matrices containing details about the new NBOs.
C                However, do not mark the 2-deg NBOs as symmetry adapted.
C
C
                 B (1,D1) = R11
                 B (2,D1) = R21
                 B (1,D2) = R12
                 B (2,D2) = R22

                 BOND1 = NBOBD (D1)
                 BOND2 = NBOBD (D2)

                 NHO1 = BDBAS (1,BOND1)
                 NHO2 = BDBAS (1,BOND2)

                 BDNBAS (BOND1) = 2
                 BDNBAS (BOND2) = 2

                 BDBAS (1,BOND1) = NHO1
                 BDBAS (2,BOND1) = NHO2
                 BDBAS (1,BOND2) = NHO1
                 BDBAS (2,BOND2) = NHO2

             END IF

         END IF
C
C
C             ...handle 3-fold degeneracy, if any.
C
C
         IF (NDEG.EQ.3) THEN

             T1 = OFFX + OFFA + FIRST
             T2 = T1 + 1
             T3 = T2 + 1

             CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                  ( NBAS,NBAS,
     +                    NBAS,NBAS,
     +                    NBAS,NDEG,
     +                    NBAS,NDEG,
     +                    NBAS,
     +                    NBAS,NBAS,3,
     +                    0,0,
     +                    SAVEC,SAVEP,SAVEC,
     +                    C,P,C (1,T1),
     +                    XVEC,
     +
     +                              XMAT )
     +
     +
             SYMSIZE = 4
             SYMSETS = NTETRA

             NSYM = SYMSIZE * SYMSETS

             DO N = 1,NSYM
                SYMCEN (N) = PLATO (N)
             END DO

             CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                  ( NBAS,NATOM,3,
     +                    NBOSIZE,
     +                    OFFLP,OFFBD,OFFEP,
     +                    OFFAB,OFFRY,OFFRR,
     +                    CENTER,
     +                    SYMSIZE,SYMSETS,SYMCEN,
     +                    SYMNBO,
     +                    QSYMACC,LSYMACC,
     +                    BDNCEN,BDCEN,
     +                    NBOBD,
     +                    SYMCRIT,
     +                    XMAT,
     +                    LOCAL,
     +                    Q,
     +
     +                             ZEROSYM,
     +                             NROT,
     +                             ROTIDX )
     +
     +
             ROT3DEG = .NOT.ZEROSYM

             IF (ROT3DEG) THEN
                 IF (NROT.EQ.0) THEN
                     WRITE (*,*) ' No 3-deg tgroup NBO rot set! '
                     WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
                     WRITE (1,*) ' No 3-deg tgroup NBO rot set! '
                     WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
                     STOP
                 END IF
             ELSE
                 WRITE (*,*) ' Could not rotate 3-deg tgroup NBOs! '
                 WRITE (*,*) ' nlo__rotate_atom_tgroup_nbo '
                 WRITE (1,*) ' Could not rotate 3-deg tgroup NBOs! '
                 WRITE (1,*) ' nlo__rotate_atom_tgroup_nbo '
                 RETURN
             END IF

             DO N = 1,NROT
                IDX = ROTIDX (N)
                XMAT (N,1) = XMAT (IDX,1)
                XMAT (N,2) = XMAT (IDX,2)
                XMAT (N,3) = XMAT (IDX,3)
             END DO

             CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                  ( 6,
     +                    ' P rotation matrix (atom,tgroup) ',
     +                    NBAS,3,
     +                    NROT,3,
     +                    XMAT )
     +
     +
C
C
C             ...find the tgroup rotation matrix.
C
C
             NORMLZE = .TRUE.
             MAXIMZE = .TRUE.

             CALL  MAT__ORTHOTRAN_MINIMAX_1_NORM
     +
     +                  ( NBAS,3,
     +                    ROTSIZE,ROTSIZE,
     +                    NBAS,NBAS,
     +                    NROT,3,
     +                    NORMLZE,
     +                    MAXIMZE,
     +                    XMAT,
     +                    XVEC (1),
     +                    XVEC (NBAS+1),
     +
     +                             NITER,
     +                             ROT )
     +
     +
             CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                  ( 6,
     +                    ' ROT matrix (atom,tgroup) ',
     +                    ROTSIZE,ROTSIZE,
     +                    3,3,
     +                    ROT )
     +
     +
C
C
C             ...form the new rotated tgroup NBO coefficients in
C                AO basis, the NBO coefficient matrix in NHO basis
C                and the matrices containing details about the new
C                NBOs. Mark the new NBOs as symmetry adapted.
C
C
             DO N = 1,NBAS
                X =   C (N,T1) * ROT (1,1)
     +              + C (N,T2) * ROT (2,1)
     +              + C (N,T3) * ROT (3,1)
                Y =   C (N,T1) * ROT (1,2)
     +              + C (N,T2) * ROT (2,2)
     +              + C (N,T3) * ROT (3,2)
                Z =   C (N,T1) * ROT (1,3)
     +              + C (N,T2) * ROT (2,3)
     +              + C (N,T3) * ROT (3,3)
                C (N,T1) = X
                C (N,T2) = Y
                C (N,T3) = Z
             END DO

             B (1,T1) = ROT (1,1)
             B (2,T1) = ROT (2,1)
             B (3,T1) = ROT (3,1)
             B (1,T2) = ROT (1,2)
             B (2,T2) = ROT (2,2)
             B (3,T2) = ROT (3,2)
             B (1,T3) = ROT (1,3)
             B (2,T3) = ROT (2,3)
             B (3,T3) = ROT (3,3)

             BOND1 = NBOBD (T1)
             BOND2 = NBOBD (T2)
             BOND3 = NBOBD (T3)

             NHO1 = BDBAS (1,BOND1)
             NHO2 = BDBAS (1,BOND2)
             NHO3 = BDBAS (1,BOND3)

             BDNBAS (BOND1) = 3
             BDNBAS (BOND2) = 3
             BDNBAS (BOND3) = 3

             BDBAS (1,BOND1) = NHO1
             BDBAS (2,BOND1) = NHO2
             BDBAS (3,BOND1) = NHO3
             BDBAS (1,BOND2) = NHO1
             BDBAS (2,BOND2) = NHO2
             BDBAS (3,BOND2) = NHO3
             BDBAS (1,BOND3) = NHO1
             BDBAS (2,BOND3) = NHO2
             BDBAS (3,BOND3) = NHO3

             SYMNBO (T1) = 2
             SYMNBO (T2) = 2
             SYMNBO (T3) = 2

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
