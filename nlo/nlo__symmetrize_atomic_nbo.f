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
         SUBROUTINE  NLO__SYMMETRIZE_ATOMIC_NBO
     +
     +                    ( NBAS,NATOM,
     +                      NBOSIZE,ROTSIZE,
     +                      MXNBA,
     +                      NXB,NXA,
     +                      OFFX,
     +                      OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR,
     +                      ATNXB,ATXIDX,ATXOFF,ATDONE,
     +                      BASBEG,BASEND,
     +                      COLMAP,SYMMAP,
     +                      IDXSYM,
     +                      INITSYM,ATOMSYM,AXISSYM,HIGHSYM,
     +                      DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                      PLATONIC,
     +                      RING,NRING,RINGSZ,
     +                      NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                      PLATO,
     +                      SYMCEN,
     +                      BDNCEN,BDCEN,
     +                      HYB,
     +                      ROT,ROTIDX,
     +                      SYMCRIT,
     +                      P,
     +                      DOMAT,
     +                      SHALF,
     +                      XVEC,
     +                      XMAT,
     +
     +                              SYMNBO,
     +                              SLTYPE,
     +                              NBOBD,
     +                              BDNBAS,
     +                              BDBAS,
     +                              LOCAL,
     +                              B,W,Q,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__SYMMETRIZE_ATOMIC_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine symmetrizes the NBOs.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    ROTSIZE      =  maximum # of NHOs that will be
C                                    allowed to form NBO hybrids for
C                                    rotationally invariant atomic NBO
C                                    construction.
C                    MXNBA        =  maximum # of NHOs per atom.
C                    NXB          =  total # of atomic x-type NBOs
C                                    (x = Core,Rydberg,Lone-pair,
C                                    Empty-pair).
C                    NXA          =  total # of atomic x-type NBO atoms.
C                    OFFX         =  absolute offset index for present
C                                    x-type NBOs. This index indicates
C                                    the totallity of NBOs before the
C                                    x-type NBO set.
C                    OFFxy        =  absolute offset indices for x-type
C                                    NBOs. These indices indicate the
C                                    totallity of x-type NBOs before the
C                                    x-type NBO set (xy:LP=Lone-pair,
C                                    BD=Bond, EP=Empty-pair, AB=Anti-
C                                    bond, RY,RR=Rydberg types)
C                    ATNXB (A)    =  # of atomic x-type NBOs on x-type
C                                    NBO atom A.
C                    ATXIDX (A)   =  atomic index for x-type NBO atom A.
C                    ATXOFF (A)   =  index offset for atomic NBOs for
C                                    x-type NBO atom A. This index is
C                                    equal to the total number of atomic
C                                    x-type NBOs on all x-type NBO atoms
C                                    preceeding x-type NBO atom A.
C                    ATDONE (A)   =  will be used as an indicator if
C                                    x-type NBO atom A has been checked
C                                    for symmetry or not.
C                    BASBEG (A)   =  first basis index number for atom A
C                    BASEND (A)   =  last basis index number for atom A
C                    COLMAP (I)   =  will contain the column map in the
C                                    active sense for reordering the
C                                    symmetrized atomic NBOs by weight
C                                    in decreasing order. COLMAP (I)
C                                    contains the position index of the
C                                    I-th symmetrized atomic NBO in the
C                                    final respective order. This array
C                                    will also be used for other
C                                    purposes as well, when rotating
C                                    NBOs on axially or atomic symmetry
C                                    related centers.
C                    SYMMAP (A,B) =  has been set equal to 1 if atoms
C                                    A and B were found to be symmetry
C                                    related. A value of 0 indicates
C                                    no symmetry relation.
C                    IDXSYM       =  array that will be used to store
C                                    symmetry related atomic indices.
C                    INITSYM      =  is true, if the routine has to
C                                    perform the initial steps, which
C                                    consists in setting the initial
C                                    values for SYMNBO (see below).
C                    ATOMSYM      =  is true, if the routine has to deal
C                                    with a single symmetry center.
C                    AXISSYM      =  is true, if the routine has to deal
C                                    with two symmetry related centers.
C                    HIGHSYM      =  is true, if the routine has to deal
C                                    with three or more symmetry related
C                                    centers.
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
C                    SYMCEN       =  array that will be used to store
C                                    symmetry related atomic indices.
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
C                    DOMAT (A,B)  =  distance order matrix containing
C                                    'distance' info between atoms A
C                                    and B. Will be used when attempting
C                                    NBO symmetrization in case of
C                                    octahedral and icosahedral groups.
C                    SHALF        =  full NBAS x NBAS square root of
C                                    the overlap matrix in AO basis.
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
C                    SLTYPE (I)   =  initial integer vector indicating
C                                    the smallest angular momentum
C                                    component present in the I-th NBO.
C                                    Useful for deciding which atomic
C                                    NBOs to use for hybrid formation.
C                                    Once an atomic NBO has been used
C                                    for hybrid formation, a -1 will
C                                    be placed at the corresponding
C                                    vector position and the NBO cannot
C                                    be used again for hybrid formation.
C                    NBOBD (I)    =  initial array that contains the
C                                    NHO bond index number for the I-th
C                                    NBO. This array is the handle for
C                                    accessing and modifying info of
C                                    the NHO bonds sitting in the arrays
C                                    BDNCEN,BDCEN,BDBAS and BDOCC.
C                    BDNBAS (J)   =  initial # of basis functions (NHOs)
C                                    for J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  initial I-th global basis (NHO)
C                                    index for J-th bond. The bond order
C                                    is NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    LOCAL        =  atomic localization content for
C                                    all NBOs before symmetrization.
C                    B (I,J)      =  initial NBOSIZE x NBAS matrix
C                                    containing the newly J-th symmetry
C                                    adapted x-type NBO expansion
C                                    coefficients in terms of the I-th
C                                    atomic NHOs forming the J-th
C                                    symmetry adapted x-type NBO.
C                    W            =  complete NBO weight vector before
C                                    symmetrization.
C                    Q            =  complete NBO occupation interaction
C                                    order vector before symmetrization.
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
C                    SLTYPE (I)   =  updated integer vector indicating
C                                    the smallest angular momentum
C                                    component present in the I-th NBO.
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
C                    NBOBD (I)    =  modified array that contains the
C                                    NHO bond index number for the I-th
C                                    NBO. This array is the handle for
C                                    accessing and modifying info of
C                                    the NHO bonds sitting in the arrays
C                                    BDNCEN,BDCEN,BDBAS and BDOCC.
C                    BDBAS (I,J)  =  modified I-th global basis (NHO)
C                                    index for J-th bond. The bond order
C                                    is NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    LOCAL        =  atomic localization content for
C                                    all NBOs after symmetrization of
C                                    the x-type NBOs treated here.
C                    B (I,J)      =  modified NBOSIZE x NBAS matrix
C                                    containing the newly J-th symmetry
C                                    adapted x-type NBO expansion
C                                    coefficients in terms of the I-th
C                                    atomic NHOs forming the J-th
C                                    symmetry adapted x-type NBO.
C                    W            =  complete NBO weight vector after
C                                    symmetrization of the x-type NBOs
C                                    treated here.
C                    Q            =  complete NBO occupation interaction
C                                    order vector after symmetrization
C                                    of the x-type NBOs treated here.
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

         LOGICAL     INITSYM,ATOMSYM,AXISSYM,HIGHSYM
         LOGICAL     C3AXIS,C4AXIS
         LOGICAL     CHECK
         LOGICAL     DEG2,DEG3,DEG4,DEG5,DEGX
         LOGICAL     LSYM,QSYM
         LOGICAL     LTRG,UTRG
         LOGICAL     PLATONIC
         LOGICAL     PROCEED
         LOGICAL     SAVEP
         LOGICAL     TGROUP,OGROUP,IGROUP
         LOGICAL     TOTAL
         LOGICAL     QUSE

         INTEGER     ATOMA,ATOMB
         INTEGER     CENTER
         INTEGER     FIRST
         INTEGER     I
         INTEGER     IDXA,IDXB
         INTEGER     MDEG
         INTEGER     MXNBA
         INTEGER     NBAS,NATOM
         INTEGER     NBO
         INTEGER     NBOSIZE,ROTSIZE
         INTEGER     NDEG
         INTEGER     NSYM
         INTEGER     NTETRA,NCUBE,NOCTA,NDODECA,NICOSA
         INTEGER     NXB,NXA
         INTEGER     NXBA
         INTEGER     OFFA,OFFX,OFFT
         INTEGER     OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR

         INTEGER     ATDONE  (1:NXA  )
         INTEGER     ATNXB   (1:NXA  )
         INTEGER     ATXIDX  (1:NXA  )
         INTEGER     ATXOFF  (1:NXA  )
         INTEGER     BASBEG  (1:NATOM)
         INTEGER     BASEND  (1:NATOM)
         INTEGER     BDNBAS  (1:NBAS )
         INTEGER     BDNCEN  (1:NBAS )
         INTEGER     COLMAP  (1:MXNBA)
         INTEGER     IDXSYM  (1:NATOM)
         INTEGER     NBOBD   (1:NBAS )
         INTEGER     NRING   (1:NATOM)
         INTEGER     RINGSZ  (1:NATOM)
         INTEGER     ROTIDX  (1:NBAS )
         INTEGER     SLTYPE  (1:NBAS )
         INTEGER     SYMCEN  (1:NATOM)
         INTEGER     SYMNBO  (1:NBAS )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS )
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS )
         INTEGER     PLATO   (1:NATOM  ,1:5    )
         INTEGER     RING    (1:NATOM  ,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM  ,1:NATOM)

         DOUBLE PRECISION  DSYMACC,LSYMACC,PSYMACC,QSYMACC
         DOUBLE PRECISION  LOLD,LNEW
         DOUBLE PRECISION  QOLD,QNEW
         DOUBLE PRECISION  SYMCRIT
         DOUBLE PRECISION  ZERO,ONE,TWO,FOUR

         DOUBLE PRECISION  LOCAL  (1:NBAS     )
         DOUBLE PRECISION  Q      (1:NBAS     )
         DOUBLE PRECISION  W      (1:NBAS     )
         DOUBLE PRECISION  XVEC   (1:NBAS+NBAS)

         DOUBLE PRECISION  B      (1:NBOSIZE,1:NBAS   )
         DOUBLE PRECISION  C      (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  HYB    (1:ROTSIZE,1:ROTSIZE)
         DOUBLE PRECISION  P      (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  ROT    (1:ROTSIZE,1:ROTSIZE)
         DOUBLE PRECISION  SHALF  (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  DOMAT  (1:NATOM  ,1:NATOM  )
         DOUBLE PRECISION  XMAT   (1:NBAS   ,1:NBAS   )

C         SAVE  TGROUP
C         SAVE  OGROUP
C         SAVE  IGROUP

         PARAMETER   (ZERO   =  0.D0 )
         PARAMETER   (ONE    =  1.D0 )
         PARAMETER   (TWO    =  2.D0 )
         PARAMETER   (FOUR   =  4.D0 )
C
C
C------------------------------------------------------------------------
C
C
C             ...perform initial step, if requested. This consists
C                of one major step:
C
C                    1) identify those x-type atomic NBOs which are
C                       non-degenerate and are thus known to be of
C                       proper symmetry already. These NBOs can be
C                       marked as usable for the symmetry reorientation
C                       procedure.
C
C
         IF (INITSYM) THEN

             DO ATOMA = 1,NXA

                IDXA = ATXIDX (ATOMA)
                OFFA = ATXOFF (ATOMA)
                NXBA = ATNXB  (ATOMA)

                OFFT = OFFX + OFFA

                NDEG = 0
                QOLD = Q     (OFFT+1)
                LOLD = LOCAL (OFFT+1)

                DO I = 1,NXBA
                   QNEW = Q     (OFFT+I)
                   LNEW = LOCAL (OFFT+I)

                   QUSE = QNEW .GT. QSYMACC
                   QSYM = DABS (QOLD - QNEW) .LT. QSYMACC
                   LSYM = DABS (LOLD - LNEW) .LT. LSYMACC

                   IF (QSYM .AND. LSYM) THEN
                       IF (QUSE) THEN
                           NDEG = NDEG + 1
                       END IF
                   ELSE
                       IF (NDEG.EQ.1) THEN
                           SYMNBO (OFFT+I-1) = 2
                       END IF
                       IF (QUSE) THEN
                           NDEG = 1
                           QOLD = QNEW
                           LOLD = LNEW
                       ELSE
                           NDEG = 0
                       END IF
                   END IF

                   IF (I.EQ.NXBA .AND. NDEG.EQ.1) THEN
                       SYMNBO (OFFT+I) = 2
                   END IF

                END DO
             END DO

             RETURN

         END IF
C
C
C             ...outer loop over all atomic centers.
C
C
         DO ATOMA = 1,NXA
            ATDONE (ATOMA) = 0
         END DO

         DO ATOMA = 1,NXA

            CHECK = ATDONE (ATOMA) .EQ. 0
C
C
C             ...find total # of symmetry related atoms corresponding
C                to current checked atom.
C
C
            IF (CHECK) THEN

                IDXA = ATXIDX (ATOMA)
                OFFA = ATXOFF (ATOMA)
                NXBA = ATNXB  (ATOMA)

                OFFT = OFFX + OFFA

                NSYM = 0
                DO ATOMB = 1,NXA
                   IDXB = ATXIDX (ATOMB)
                   IF (SYMMAP (IDXB,IDXA) .EQ. 1) THEN
                       NSYM = NSYM + 1
                       IDXSYM (NSYM) = ATOMB
                       ATDONE (ATOMB) = 1
                   END IF
                END DO
C
C
C             ...proceed with the high symmetry case, which deals
C                with the NSYM > 2 cases. 2-fold degenerate levels
C                will be identified by both identical interaction
C                orders and atomic localization content.
C
C
                IF (HIGHSYM .AND. NSYM.GT.2) THEN

                    NDEG = 0
                    QOLD = Q     (OFFT+1)
                    LOLD = LOCAL (OFFT+1)

                    DO I = 1,NXBA
                       QNEW = Q     (OFFT+I)
                       LNEW = LOCAL (OFFT+I)

                       QUSE = QNEW .GT. QSYMACC
                       QSYM = DABS (QOLD - QNEW) .LT. QSYMACC
                       LSYM = DABS (LOLD - LNEW) .LT. LSYMACC

                       IF (QSYM .AND. LSYM) THEN
                           IF (QUSE) THEN
                               NDEG = NDEG + 1
                           END IF
                       ELSE
                           IF (QUSE) THEN
                               NDEG = 1
                               QOLD = QNEW
                               LOLD = LNEW
                           ELSE
                               NDEG = 0
                           END IF
                       END IF

                       IF (NDEG.GT.2) THEN
                           WRITE (*,*) ' Degeneracy > 2-fold! '
                           WRITE (*,*) ' Atom index = ',IDXA
                           WRITE (*,*) ' Order,Locality = ',QOLD,LOLD
                           WRITE (*,*) ' nlo__symmetrize_atomic_nbo '
                           WRITE (1,*) ' Degeneracy > 2-fold! '
                           WRITE (1,*) ' Atom index = ',IDXA
                           WRITE (1,*) ' Order,Locality = ',QOLD,LOLD
                           WRITE (1,*) ' nlo__symmetrize_atomic_nbo '
                           STOP
                       END IF

                       IF (NDEG.EQ.2) THEN

                           FIRST = I - 1

                           WRITE (*,*) ' HIGH SYMMETRY case! '
                           WRITE (*,*) ' Weight = ',W (OFFT+FIRST)
                           WRITE (*,*) ' Order = ',QOLD
                           WRITE (*,*) ' Loc = ',LOLD
                           WRITE (*,*) ' 1st function = ',FIRST
                           WRITE (*,*) ' High NSYM = ',NSYM

                           IF (NSYM.EQ.4) THEN

                               CALL  NLO__ROTATE_TGROUP_NBO
     +
     +                                    ( NBAS,
     +                                      NBOSIZE,
     +                                      NXB,NXA,NXBA,
     +                                      ATXOFF,
     +                                      LSYMACC,QSYMACC,
     +                                      FIRST,
     +                                      IDXSYM,
     +                                      NBOBD (OFFX+1),
     +                                      P,
     +                                      LOCAL (OFFT+1),
     +                                      W (OFFT+1),
     +                                      Q (OFFT+1),
     +                                      XVEC,
     +                                      XMAT,
     +
     +                                             SYMNBO (OFFX+1),
     +                                             SLTYPE (OFFX+1),
     +                                             BDNBAS,
     +                                             BDBAS,
     +                                             B (1,OFFX+1),
     +                                             C (1,OFFX+1) )
     +
     +
                               TGROUP = .TRUE.

                           ELSE IF (NSYM.EQ.6 .OR. NSYM.EQ.8) THEN

                               C3AXIS = NSYM.EQ.8
                               C4AXIS = NSYM.EQ.6

                               CALL  NLO__ROTATE_OGROUP_NBO
     +
     +                                    ( NBAS,NATOM,NSYM,
     +                                      NBOSIZE,
     +                                      NXB,NXA,NXBA,
     +                                      ATXIDX,ATXOFF,
     +                                      DSYMACC,LSYMACC,QSYMACC,
     +                                      FIRST,
     +                                      IDXSYM,
     +                                      C3AXIS,C4AXIS,
     +                                      NBOBD (OFFX+1),
     +                                      P,
     +                                      LOCAL (OFFT+1),
     +                                      W (OFFT+1),
     +                                      Q (OFFT+1),
     +                                      DOMAT,
     +                                      XVEC,
     +                                      XMAT,
     +
     +                                             SYMNBO (OFFX+1),
     +                                             SLTYPE (OFFX+1),
     +                                             BDNBAS,
     +                                             BDBAS,
     +                                             B (1,OFFX+1),
     +                                             C (1,OFFX+1) )
     +
     +
                               OGROUP = .TRUE.

                           ELSE IF (NSYM.EQ.12) THEN

                               WRITE (*,*) ' I,Ih-group present! '
                               WRITE (*,*) ' Not implemented yet! '
                               WRITE (*,*) ' Symmetrization skipped... '
                               WRITE (1,*) ' I,Ih-group present! '
                               WRITE (1,*) ' Not implemented yet! '
                               WRITE (1,*) ' Symmetrization skipped... '

                               IGROUP = .TRUE.

                           ELSE IF (NSYM.EQ.20) THEN

                               WRITE (*,*) ' I,Ih-group present! '
                               WRITE (*,*) ' Not implemented yet! '
                               WRITE (*,*) ' Symmetrization skipped... '
                               WRITE (1,*) ' I,Ih-group present! '
                               WRITE (1,*) ' Not implemented yet! '
                               WRITE (1,*) ' Symmetrization skipped... '

                               IGROUP = .TRUE.

                           ELSE

                               WRITE (*,*) ' # of symcen ne 4,6,8,12,20'
                               WRITE (*,*) ' Accidental degeneracy ? '
                               WRITE (*,*) ' nlo__symmetrize_atomic_nbo'
                               WRITE (1,*) ' # of symcen ne 4,6,8,12,20'
                               WRITE (1,*) ' Accidental degeneracy ? '
                               WRITE (1,*) ' nlo__symmetrize_atomic_nbo'
                               STOP

                           END IF

                       END IF
                    END DO

                END IF
C
C
C             ...proceed with the axial symmetry case, which deals
C                with the NSYM = 2 cases. The following point groups
C                are possible here:
C
C                        Ci,Cs,C2,C2v,Cnh,Dn,Dnh,Dnd,Sn
C
C                Of these, only the following will lead to 2-fold
C                degenerate levels (i.e. only those groups will allow
C                for site-groups of Cn,Cnv,n>=3):
C
C                            Cnh,Dn,Dnh,Dnd  n>=3  and Sn  n>=6
C
C                2-fold degenerate levels will be identified by both
C                identical interaction orders and atomic localization
C                content. Note also that we cannot have the axial
C                symmetry case, if a high symmetry group was detected
C                before.
C
C
                IF (AXISSYM .AND. NSYM.EQ.2) THEN

                    IF (PLATONIC) THEN
                        WRITE (*,*) ' Axial symmetry conflict! '
                        WRITE (*,*) ' PLATONIC = ',PLATONIC
                        WRITE (*,*) ' nlo__symmetrize_atomic_nbo'
                        WRITE (1,*) ' Axial symmetry conflict! '
                        WRITE (1,*) ' PLATONIC = ',PLATONIC
                        WRITE (1,*) ' nlo__symmetrize_atomic_nbo'
                        STOP
                    END IF

                    NDEG = 0
                    QOLD = Q     (OFFT+1)
                    LOLD = LOCAL (OFFT+1)

                    DO I = 1,NXBA
                       QNEW = Q     (OFFT+I)
                       LNEW = LOCAL (OFFT+I)

                       QUSE = QNEW .GT. QSYMACC
                       QSYM = DABS (QOLD - QNEW) .LT. QSYMACC
                       LSYM = DABS (LOLD - LNEW) .LT. LSYMACC

                       IF (QSYM .AND. LSYM) THEN
                           IF (QUSE) THEN
                               NDEG = NDEG + 1
                           END IF
                       ELSE
                           IF (QUSE) THEN
                               NDEG = 1
                               QOLD = QNEW
                               LOLD = LNEW
                           ELSE
                               NDEG = 0
                           END IF
                       END IF

                       IF (NDEG.GT.2) THEN
                           WRITE (*,*) ' Degeneracy > 2-fold! '
                           WRITE (*,*) ' Atom index = ',IDXA
                           WRITE (*,*) ' Order,Locality = ',QOLD,LOLD
                           WRITE (*,*) ' nlo__symmetrize_atomic_nbo '
                           WRITE (1,*) ' Degeneracy > 2-fold! '
                           WRITE (1,*) ' Atom index = ',IDXA
                           WRITE (1,*) ' Order,Locality = ',QOLD,LOLD
                           WRITE (1,*) ' nlo__symmetrize_atomic_nbo '
                           STOP
                       END IF

                       IF (NDEG.EQ.2) THEN

                           FIRST = I - 1
                           PROCEED = SYMNBO (OFFT+FIRST) .NE. 2

                           IF (PROCEED) THEN
                               WRITE (*,*) ' AXIAL case! '
                               WRITE (*,*) ' Weight = ',W (OFFT+FIRST)
                               WRITE (*,*) ' Order = ',Q (OFFT+FIRST)
                               WRITE (*,*) ' Loc = ',LOCAL (OFFT+FIRST)
                               WRITE (*,*) ' Degeneracy = ',NDEG
                               WRITE (*,*) ' 1st function = ',FIRST
                               WRITE (*,*) ' Center # 1 = ',IDXSYM (1)
                               WRITE (*,*) ' Center # 2 = ',IDXSYM (2)

                               CALL  NLO__ROTATE_AXIAL_NBO
     +
     +                                    ( NBAS,NATOM,
     +                                      NBOSIZE,ROTSIZE,
     +                                      NXA,NXBA,
     +                                      OFFX,
     +                                      OFFLP,OFFBD,OFFEP,
     +                                      OFFAB,OFFRY,OFFRR,
     +                                      ATXIDX,ATXOFF,
     +                                      IDXSYM,
     +                                      RING,NRING,RINGSZ,
     +                                      SYMCEN,
     +                                      LSYMACC,QSYMACC,
     +                                      FIRST,
     +                                      COLMAP,
     +                                      NBOBD,
     +                                      BDNCEN,BDCEN,
     +                                      HYB,
     +                                      ROT,ROTIDX,
     +                                      SYMCRIT,
     +                                      P,
     +                                      LOCAL,
     +                                      W,Q,
     +                                      XVEC,
     +                                      XMAT,
     +
     +                                              SYMNBO,
     +                                              SLTYPE,
     +                                              BDNBAS,
     +                                              BDBAS,
     +                                              B,
     +                                              C )
     +
     +
                           END IF
                       END IF
                    END DO

                END IF
C
C
C             ...proceed with the atomic symmetry case, which deals
C                with the NSYM = 1 cases. 2- or higher-fold degenerate
C                levels will be identified by both identical interaction
C                orders and atomic localization content. If 2- or
C                higher-fold degenerate levels are found, this means
C                that either the atomic center is located at the center
C                point of the point group or it is located on a C axis
C                of order =< 3. Note that if at least one 3- or higher-
C                fold degenerate level is found, then the atomic center
C                must be located at the center point of the point group
C                and the group can only be of T-,O- or I-type.
C
C
                IF (ATOMSYM .AND. NSYM.EQ.1) THEN

                    CENTER = ATXIDX (IDXSYM (1))

                    DEG2 = .FALSE.
                    DEG3 = .FALSE.
                    DEG4 = .FALSE.
                    DEG5 = .FALSE.
                    DEGX = .FALSE.
                    NDEG = 0
                    QOLD = Q     (OFFT+1)
                    LOLD = LOCAL (OFFT+1)

                    DO I = 1,NXBA
                       QNEW = Q     (OFFT+I)
                       LNEW = LOCAL (OFFT+I)

                       QUSE = QNEW .GT. QSYMACC
                       QSYM = DABS (QOLD - QNEW) .LT. QSYMACC
                       LSYM = DABS (LOLD - LNEW) .LT. LSYMACC

                       IF (QSYM .AND. LSYM) THEN
                           IF (QUSE) THEN
                               NDEG = NDEG + 1
                               IF (I.EQ.NXBA) THEN
                                   DEG2 = NDEG.EQ.2
                                   DEG3 = NDEG.EQ.3
                                   DEG4 = NDEG.EQ.4
                                   DEG5 = NDEG.EQ.5
                                   DEGX = NDEG.GT.5
                                   MDEG = NDEG
                                   FIRST = I - MDEG + 1
                               END IF
                           END IF
                       ELSE
                           DEG2 = NDEG.EQ.2
                           DEG3 = NDEG.EQ.3
                           DEG4 = NDEG.EQ.4
                           DEG5 = NDEG.EQ.5
                           DEGX = NDEG.GT.5
                           MDEG = NDEG
                           FIRST = I - MDEG
                           IF (QUSE) THEN
                               NDEG = 1
                               QOLD = QNEW
                               LOLD = LNEW
                           ELSE
                               NDEG = 0
                           END IF
                       END IF

                       PROCEED =  (DEG2 .OR. DEG3 .OR. DEG4 .OR. DEG5)
     +                                           .AND.
     +                               (SYMNBO (OFFT+FIRST) .NE. 2)

                       IF (PROCEED) THEN
                           WRITE (*,*) ' ATOMIC case! '
                           WRITE (*,*) ' Weight = ',W (OFFT+FIRST)
                           WRITE (*,*) ' Order = ',Q (OFFT+FIRST)
                           WRITE (*,*) ' Loc = ',LOCAL (OFFT+FIRST)
                           WRITE (*,*) ' Degeneracy = ',MDEG
                           WRITE (*,*) ' 1st function = ',FIRST

                           IF (NTETRA.GT.0) THEN

                               CALL  NLO__ROTATE_ATOM_TGROUP_NBO
     +
     +                                    ( NBAS,NATOM,MDEG,
     +                                      NBOSIZE,ROTSIZE,
     +                                      NXBA,
     +                                      OFFX,OFFA,
     +                                      OFFLP,OFFBD,OFFEP,
     +                                      OFFAB,OFFRY,OFFRR,
     +                                      CENTER,
     +                                      NTETRA,
     +                                      PLATO (1,1),
     +                                      SYMCEN,
     +                                      LSYMACC,QSYMACC,
     +                                      FIRST,
     +                                      COLMAP,
     +                                      NBOBD,
     +                                      BDNCEN,BDCEN,
     +                                      HYB,
     +                                      ROT,ROTIDX,
     +                                      SYMCRIT,
     +                                      P,
     +                                      LOCAL,
     +                                      W,Q,
     +                                      XVEC,
     +                                      XMAT,
     +
     +                                               SYMNBO,
     +                                               SLTYPE,
     +                                               BDNBAS,
     +                                               BDBAS,
     +                                               B,
     +                                               C )
     +
     +
                           ELSE IF ((NCUBE.GT.0).OR.(NOCTA.GT.0)) THEN

                               CALL  NLO__ROTATE_ATOM_OGROUP_NBO
     +
     +                                    ( NBAS,NATOM,MDEG,
     +                                      NBOSIZE,ROTSIZE,
     +                                      NXBA,
     +                                      OFFX,OFFA,
     +                                      OFFLP,OFFBD,OFFEP,
     +                                      OFFAB,OFFRY,OFFRR,
     +                                      CENTER,
     +                                      NCUBE,NOCTA,
     +                                      PLATO (1,2),
     +                                      SYMCEN,
     +                                      LSYMACC,QSYMACC,
     +                                      FIRST,
     +                                      COLMAP,
     +                                      NBOBD,
     +                                      BDNCEN,BDCEN,
     +                                      HYB,
     +                                      ROT,ROTIDX,
     +                                      SYMCRIT,
     +                                      P,
     +                                      LOCAL,
     +                                      W,Q,
     +                                      XVEC,
     +                                      XMAT,
     +
     +                                               SYMNBO,
     +                                               SLTYPE,
     +                                               BDNBAS,
     +                                               BDBAS,
     +                                               B,
     +                                               C )
     +
     +
                           ELSE IF (NDODECA.GT.0.OR.NICOSA.GT.0) THEN

                               WRITE (*,*) ' I,Ih-group present! '
                               WRITE (*,*) ' Not implemented yet! '
                               WRITE (*,*) ' Symmetrization skipped... '
                               WRITE (1,*) ' I,Ih-group present! '
                               WRITE (1,*) ' Not implemented yet! '
                               WRITE (1,*) ' Symmetrization skipped... '

                           ELSE IF (PLATONIC) THEN

                               WRITE (*,*) ' Platonic sym present! '
                               WRITE (*,*) ' But no surroundings! '
                               WRITE (*,*) ' Symmetrization skipped... '
                               WRITE (1,*) ' Platonic sym present! '
                               WRITE (1,*) ' But no surroundings! '
                               WRITE (1,*) ' Symmetrization skipped... '

                           ELSE

                               CALL  NLO__ROTATE_ATOM_NBO
     +
     +                                    ( NBAS,NATOM,MDEG,
     +                                      NBOSIZE,ROTSIZE,
     +                                      NXBA,
     +                                      OFFX,OFFA,
     +                                      OFFLP,OFFBD,OFFEP,
     +                                      OFFAB,OFFRY,OFFRR,
     +                                      CENTER,
     +                                      RING,NRING,RINGSZ,
     +                                      SYMCEN,
     +                                      LSYMACC,QSYMACC,
     +                                      FIRST,
     +                                      COLMAP,
     +                                      NBOBD,
     +                                      BDNCEN,BDCEN,
     +                                      HYB,
     +                                      ROT,ROTIDX,
     +                                      SYMCRIT,
     +                                      P,
     +                                      LOCAL,
     +                                      W,Q,
     +                                      XVEC,
     +                                      XMAT,
     +
     +                                              SYMNBO,
     +                                              SLTYPE,
     +                                              BDNBAS,
     +                                              BDBAS,
     +                                              B,
     +                                              C )
     +
     +
                           END IF

                           DEG2 = .FALSE.
                           DEG3 = .FALSE.
                           DEG4 = .FALSE.
                           DEG5 = .FALSE.

                       END IF

                    END DO

                END IF
C
C
C             ...check next atom.
C
C
            END IF
         END DO
C
C
C             ...establish the final atomic NBO weights, localization
C                contents and interaction orders, and reorder the atomic
C                NBOs according to decreasing weight.
C
C
         LTRG = .FALSE.
         UTRG = .FALSE.
         SAVEP = .TRUE.
         TOTAL = .TRUE.

         DO ATOMA = 1,NXA

            IDXA = ATXIDX (ATOMA)
            OFFA = ATXOFF (ATOMA)
            NXBA = ATNXB  (ATOMA)

            OFFT = OFFX + OFFA

            CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                 ( NBAS,NXBA,
     +                   NBAS,NBAS,
     +                   NBAS,NXBA,
     +                   NXBA,
     +                   NXBA,NBAS,
     +                   0,
     +                   SAVEP,LTRG,UTRG,
     +                   C (1,OFFT+1),
     +                   P,
     +
     +                           W (OFFT+1),
     +
     +                   XMAT )
     +
     +
            DO NBO = OFFT+1,OFFT+NXBA

               CALL  NLO__EXTRACT_ATOMIC_CONTENT
     +
     +                    ( NBAS,NATOM,
     +                      1,
     +                      IDXA,
     +                      BASBEG,BASEND,
     +                      TOTAL,
     +                      C (1,NBO),
     +                      SHALF,
     +
     +                               LOCAL (NBO) )
     +
     +
               CALL  NLO__CALC_INTERACTION_ORDER
     +
     +                    ( NBAS,NBAS,
     +                      NBAS,NBAS,
     +                      NBAS,
     +                      NBAS,NBAS,
     +                      NBAS,
     +                      NBO,
     +                      C,P,C (1,NBO),
     +                      XVEC (1),XVEC (NBAS+1),
     +
     +                               Q (NBO) )
     +
     +
            END DO
C
C
C             ...reorder curent atomic NBOs according to decreasing
C                weight. Reorder also all the data associated with
C                these NBOs. Obsolete array ROTIDX will be used
C                as an integer scratch vector.
C
C
            CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +                 ( NXBA,NXBA,
     +                   1,NXBA,
     +                   .FALSE.,.FALSE.,
     +                   0,
     +                   W (OFFT+1),
     +
     +                           ROTIDX )
     +
     +
            DO I = 1,NXBA
               COLMAP (ROTIDX (I)) = I
            END DO

            CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +                 ( NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   COLMAP,
     +                   XVEC,
     +
     +                           LOCAL (OFFT+1) )
     +
     +
            CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +                 ( NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   COLMAP,
     +                   XVEC,
     +
     +                           Q (OFFT+1) )
     +
     +
            CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +                 ( NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   COLMAP,
     +                   XVEC,
     +
     +                           W (OFFT+1) )
     +
     +
            CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +                 ( NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   COLMAP,
     +                   ROTIDX,
     +
     +                           NBOBD (OFFT+1) )
     +
     +
            CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +                 ( NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   COLMAP,
     +                   ROTIDX,
     +
     +                           SYMNBO (OFFT+1) )
     +
     +
            CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +                 ( NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   COLMAP,
     +                   ROTIDX,
     +
     +                           SLTYPE (OFFT+1) )
     +
     +
            CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +                 ( NBOSIZE,NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NBOSIZE,
     +                   NBOSIZE,NXBA,
     +                   COLMAP,
     +                   ROTIDX,
     +                   XVEC,
     +
     +                           B (1,OFFT+1) )
     +
     +
            CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +                 ( NBAS,NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   NBAS,
     +                   NBAS,NXBA,
     +                   COLMAP,
     +                   ROTIDX,
     +                   XVEC,
     +
     +                           C (1,OFFT+1) )
     +
     +
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
