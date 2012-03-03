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
         SUBROUTINE  NLO__ROTATE_ATOM_OGROUP_NBO
     +
     +                    ( NBAS,NATOM,NDEG,
     +                      NBOSIZE,ROTSIZE,
     +                      NXBA,
     +                      OFFX,OFFA,
     +                      OFFLP,OFFBD,OFFEP,
     +                      OFFAB,OFFRY,OFFRR,
     +                      CENTER,
     +                      NCUBE,NOCTA,
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
C  OPERATION   : NLO__ROTATE_ATOM_OGROUP_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine rotates a degenerate set of NBOs sitting
C                at the center of octahedral symmetry. The possible
C                background atomic sets used for performing this
C                rotation can be either chosen from the octahedrally or
C                cubically related centers, where the former posses
C                C4v and the latter C3v site symmetry. The following
C                Oh -> C4v and C3v irrep decomposition table shows what
C                will be expected:
C
C
C                         Oh   |    C4v   |    C3v
C                      ------------------------------
C                         A1g  |     A1   |    A1
C                         A2g  |     B1   |    A2
C                          Eg  |  A1 + B1 |     E
C                         T1g  |  A2 + E  |  A2 + E
C                         T2g  |  B2 + E  |  A1 + E
C                         A1u  |     A2   |    A2
C                         A2u  |     B2   |    A1
C                          Eu  |  A2 + B2 |     E
C                         T1u  |  A1 + E  |  A1 + E
C                         T2u  |  B1 + E  |  A2 + E
C
C
C                From this table we can construct two tables, showing
C                which octahedral degenerate NBO sets can have (+)
C                or cannot have (-) interactions with the octahedral
C                (C4v site symmetry) and cubical (C3v site symmetry)
C                background:
C
C
C
C                    C4v background | Oh-irrep: Eg T1g T2g Eu T1u T2u
C                  ----------------------------------------------------
C                    A1 (s,pz)      |           +   -   -  -   +   -
C                    A2 (g,h,i,...) |           -   +   -  +   -   -
C                    B1 (dx^2-y^2)  |           +   -   -  -   -   +
C                    B2 (dxy)       |           -   -   +  +   -   -
C                     E (px,py)     |           -   +   +  -   +   +
C
C
C
C                    C3v background | Oh-irrep: Eg T1g T2g Eu T1u T2u
C                  ----------------------------------------------------
C                    A1 (s,pz)      |           -   -   +  -   +   -
C                    A2 (f,g,h,...) |           -   +   -  -   -   +
C                  A1+E (sp2-hybrid)|           +   +   +  +   +   +
C
C
C
C                Both tables also show what type of C4v and C3v irreps
C                will occur for specific atomic functions used on those
C                sites. For example, when s,p,d,f functions are used
C                on one C4v site we expect all C4v irreps to occur,
C                except the A2 irrep, which will only occur if g- or
C                higher-angular functions are used. Similarly A2 irreps
C                of C3v will only occur when including f- or higher
C                atomic functions.
C
C                Another table is also useful, showing decomposition
C                of atomic angular momentum levels into Oh irreps:
C
C
C                        atomic level  |          Oh irreps
C                      -----------------------------------------------
C                             s        |             A1g
C                             p        |             T1u
C                             d        |           Eg + T2g
C                             f        |        A2u + T1u + T2u
C                             g        |     A1g + Eg + T1g + T2g
C                             h        |        Eu + 2T1u + T2u
C                             i        |  A1g + A2g + Eg + T1g + 2T2g
C
C
C                We have next to distinguish two cases:
C
C                  i) 2-fold degenerate NBOs
C
C                     These can be only of Eg and Eu type and can
C                     occur if the basis set on the center atom has
C                     at least respectively d-type and h-type functions.
C                     In the majority of the cases we will thus deal
C                     with Eg type degenerate NBOs. Eu type NBOs can
C                     be checked by the odd smallest angular momentum
C                     component presence in these NBOs.
C
C                     Procedure for the Eg type NBOs:
C
C                       a) Rotate the Eg pair to one of the background
C                          sites, such that the first Eg NBO interaction
C                          with this site is either minimized (C4v
C                          octahedral A1 or B1 site symmetry) or
C                          maximized (C3v cubical A1+E sp2 hybrid site
C                          symmetry).
C
C                       b) Once appropriately rotated, form sd2
C                          hybrids to achieve octahedral symmetry.
C                           
C
C                 ii) 3-fold degenerate NBOs
C
C                     These can be of T1g,T2g,T1u and T2u type. Both
C                     background types C4v and C3v are checked for
C                     rotation, with the octahedral C4v preferred over
C                     the cubical C3v (less background points for the
C                     former!). The 1-norm of the interaction occupation
C                     matrix between the background and the 3-fold
C                     degenerate NBOs is then:
C
C                          1-norm minimized -> for C4v background
C                          1-norm maximized -> for C3v background
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
C                    NCUBE        =  contains the # of cubically related
C                                    center sets (size 8)
C                    NOCTA        =  contains the # of octahedrally
C                                    related center sets (size 6)
C                    PLATO (A,P)  =  contains all atomic indices A for
C                                    all cubically arranged centers
C                                    in groups of 8 (P = 1) and all
C                                    octahedrally arranged centers in
C                                    groups of 6 (P = 2).
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

         LOGICAL     CUBICAL
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
         INTEGER     NCUBE,NOCTA
         INTEGER     NHO1,NHO2,NHO3
         INTEGER     NDEG
         INTEGER     NITER
         INTEGER     NLO__ANGULAR_MOMENTUM_NBO
         INTEGER     NSYM
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
         INTEGER     ROTIDX  (1:NBAS )
         INTEGER     SLTYPE  (1:NBAS )
         INTEGER     SYMCEN  (1:NATOM)
         INTEGER     SYMNBO  (1:NBAS )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS )
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS )
         INTEGER     PLATO   (1:NATOM  ,1:2    )

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
C                with the octahedral NBOs, we know that deg = 2 or 3
C                always. Otherwise stop with message.
C
C
         IF (NDEG.NE.2 .AND. NDEG.NE.3) THEN
             WRITE (*,*) ' Degeneracy conflict for atomic ogroup case! '
             WRITE (*,*) ' NDEG = ',NDEG
             WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
             WRITE (1,*) ' Degeneracy conflict for atomic ogroup case! '
             WRITE (1,*) ' NDEG = ',NDEG
             WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
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
             WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
             WRITE (1,*) ' Degeneracy size too large! '
             WRITE (1,*) ' NDEG,ROTSIZE = ',NDEG,ROTSIZE
             WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
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
C             ...handle 2-fold degeneracy, if any. Try to use first
C                the octahedrally arranged centers, if any, followed
C                by the cubically arranged ones.
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
             ROT2DEG = .FALSE.

C             GOTO 1234

             IF (NOCTA.NE.0) THEN

                 SYMSIZE = 6
                 SYMSETS = NOCTA

                 NSYM = SYMSIZE * SYMSETS

                 DO N = 1,NSYM
                    SYMCEN (N) = PLATO (N,2)
                 END DO

                 CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                      ( NBAS,NATOM,2,
     +                        NBOSIZE,
     +                        OFFLP,OFFBD,OFFEP,
     +                        OFFAB,OFFRY,OFFRR,
     +                        CENTER,
     +                        SYMSIZE,SYMSETS,SYMCEN,
     +                        SYMNBO,
     +                        QSYMACC,LSYMACC,
     +                        BDNCEN,BDCEN,
     +                        NBOBD,
     +                        SYMCRIT,
     +                        XMAT,
     +                        LOCAL,
     +                        Q,
     +
     +                                 ZEROSYM,
     +                                 NROT,
     +                                 ROTIDX )
     +
     +
                 ROT2DEG = .NOT.ZEROSYM

                 IF (ROT2DEG) THEN
                     WRITE (*,*) ' Using 2-deg octahedral background! '
                     IF (NROT.EQ.0) THEN
                         WRITE (*,*) ' No 2-deg ogroup NBO rot set! '
                         WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                         WRITE (1,*) ' No 2-deg ogroup NBO rot set! '
                         WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
                         STOP
                     END IF
                 END IF
             END IF

C             GOTO 12345

 1234        IF (.NOT.ROT2DEG .AND. (NCUBE.NE.0)) THEN

                 SYMSIZE = 8
                 SYMSETS = NCUBE

                 NSYM = SYMSIZE * SYMSETS

                 DO N = 1,NSYM
                    SYMCEN (N) = PLATO (N,1)
                 END DO

                 CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                      ( NBAS,NATOM,2,
     +                        NBOSIZE,
     +                        OFFLP,OFFBD,OFFEP,
     +                        OFFAB,OFFRY,OFFRR,
     +                        CENTER,
     +                        SYMSIZE,SYMSETS,SYMCEN,
     +                        SYMNBO,
     +                        QSYMACC,LSYMACC,
     +                        BDNCEN,BDCEN,
     +                        NBOBD,
     +                        SYMCRIT,
     +                        XMAT,
     +                        LOCAL,
     +                        Q,
     +
     +                                 ZEROSYM,
     +                                 NROT,
     +                                 ROTIDX )
     +
     +
                 ROT2DEG = .NOT.ZEROSYM

                 IF (ROT2DEG) THEN
                     WRITE (*,*) ' Using 2-deg cubical background! '
                     IF (NROT.EQ.0) THEN
                         WRITE (*,*) ' No 2-deg ogroup NBO rot set! '
                         WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                         WRITE (1,*) ' No 2-deg ogroup NBO rot set! '
                         WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
                         STOP
                     END IF
                 END IF
             END IF

12345        IF (.NOT.ROT2DEG) THEN
                 WRITE (*,*) ' Could not rotate 2-deg ogroup NBOs! '
                 WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                 WRITE (1,*) ' Could not rotate 2-deg ogroup NBOs! '
                 WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
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
C                The same logic applies to the Eu-irreps, which
C                appear for the first time with h-functions and are
C                of the same form as those for the d-functions but
C                multiplied both be a x*y*z factor.
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
                     WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                     WRITE (1,*) ' Cannot form sd2 hybrids! '
                     WRITE (1,*) ' ROTSIZE = ',ROTSIZE
                     WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
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
     +                        ' P for sd2 orientation (atom,ogroup) ',
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
                 WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                 WRITE (*,*) ' O-group symmetry reduction will occur! '
                 WRITE (1,*) ' No s function for sd2 hybrid! '
                 WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
                 WRITE (1,*) ' O-group symmetry reduction will occur! '

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
C             ...handle 3-fold degeneracy, if any. Try to use first
C                the octahedrally arranged centers, if any, followed
C                by the cubically arranged ones.
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
             CUBICAL = .FALSE.
             ROT3DEG = .FALSE.

C             GOTO 2345

             IF (NOCTA.NE.0) THEN

                 SYMSIZE = 6
                 SYMSETS = NOCTA

                 NSYM = SYMSIZE * SYMSETS

                 DO N = 1,NSYM
                    SYMCEN (N) = PLATO (N,2)
                 END DO

                 CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                      ( NBAS,NATOM,3,
     +                        NBOSIZE,
     +                        OFFLP,OFFBD,OFFEP,
     +                        OFFAB,OFFRY,OFFRR,
     +                        CENTER,
     +                        SYMSIZE,SYMSETS,SYMCEN,
     +                        SYMNBO,
     +                        QSYMACC,LSYMACC,
     +                        BDNCEN,BDCEN,
     +                        NBOBD,
     +                        SYMCRIT,
     +                        XMAT,
     +                        LOCAL,
     +                        Q,
     +
     +                                 ZEROSYM,
     +                                 NROT,
     +                                 ROTIDX )
     +
     +
                 ROT3DEG = .NOT.ZEROSYM

                 IF (ROT3DEG) THEN
                     WRITE (*,*) ' Using 3-deg octahedral background! '
                     IF (NROT.EQ.0) THEN
                         WRITE (*,*) ' No 3-deg ogroup NBO rot set! '
                         WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                         WRITE (1,*) ' No 3-deg ogroup NBO rot set! '
                         WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
                         STOP
                     END IF
                 END IF
             END IF

C             GOTO 23456

 2345        IF (.NOT.ROT3DEG .AND. (NCUBE.NE.0)) THEN

                 SYMSIZE = 8
                 SYMSETS = NCUBE

                 NSYM = SYMSIZE * SYMSETS

                 DO N = 1,NSYM
                    SYMCEN (N) = PLATO (N,1)
                 END DO

                 CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                      ( NBAS,NATOM,3,
     +                        NBOSIZE,
     +                        OFFLP,OFFBD,OFFEP,
     +                        OFFAB,OFFRY,OFFRR,
     +                        CENTER,
     +                        SYMSIZE,SYMSETS,SYMCEN,
     +                        SYMNBO,
     +                        QSYMACC,LSYMACC,
     +                        BDNCEN,BDCEN,
     +                        NBOBD,
     +                        SYMCRIT,
     +                        XMAT,
     +                        LOCAL,
     +                        Q,
     +
     +                                 ZEROSYM,
     +                                 NROT,
     +                                 ROTIDX )
     +
     +
                 ROT3DEG = .NOT.ZEROSYM

                 IF (ROT3DEG) THEN
                     WRITE (*,*) ' Using 3-deg cubical background! '

                     CUBICAL = .TRUE.

                     IF (NROT.EQ.0) THEN
                         WRITE (*,*) ' No 3-deg ogroup NBO rot set! '
                         WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                         WRITE (1,*) ' No 3-deg ogroup NBO rot set! '
                         WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
                         STOP
                     END IF
                 END IF
             END IF

23456        IF (.NOT.ROT3DEG) THEN
                 WRITE (*,*) ' Could not rotate 3-deg ogroup NBOs! '
                 WRITE (*,*) ' nlo__rotate_atom_ogroup_nbo '
                 WRITE (1,*) ' Could not rotate 3-deg ogroup NBOs! '
                 WRITE (1,*) ' nlo__rotate_atom_ogroup_nbo '
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
     +                    ' P rotation matrix (atom,ogroup) ',
     +                    NBAS,3,
     +                    NROT,3,
     +                    XMAT )
     +
     +
C
C
C             ...find the ogroup rotation matrix.
C
C
             NORMLZE = .TRUE.
             MAXIMZE = CUBICAL

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
     +                    ' ROT matrix (atom,ogroup) ',
     +                    ROTSIZE,ROTSIZE,
     +                    3,3,
     +                    ROT )
     +
     +
C
C
C             ...form the new rotated ogroup NBO coefficients in
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
