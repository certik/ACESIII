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
         SUBROUTINE  NLO__ROTATE_AXIAL_NBO
     +
     +                    ( NBAS,NATOM,
     +                      NBOSIZE,ROTSIZE,
     +                      NXA,NXBA,
     +                      OFFX,
     +                      OFFLP,OFFBD,OFFEP,
     +                      OFFAB,OFFRY,OFFRR,
     +                      ATXIDX,ATXOFF,
     +                      IDXSYM,
     +                      RING,NRING,RINGSZ,
     +                      SYMCEN,
     +                      LSYMACC,QSYMACC,
     +                      NPAIR,
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
C  OPERATION   : NLO__ROTATE_AXIAL_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine rotates an axial set (two) of two-fold
C                interaction order degenerate atomic NBOs, such that
C                they become equivalent under all axial operations.
C
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
C                    NXA          =  total # of atomic x-type NBO atoms.
C                    NXBA         =  # of atomic x-type NBOs on one of
C                                    the axially related atoms.
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
C                    ATXIDX (A)   =  atomic index for x-type NBO atom A.
C                    ATXOFF (A)   =  index offset for atomic NBOs for
C                                    x-type NBO atom A. This index is
C                                    equal to the total number of atomic
C                                    x-type NBOs on all x-type NBO atoms
C                                    preceeding x-type NBO atom A.
C                    IDXSYM (1,2) =  both x-type NBO atomic indices
C                                    which are axially symmetry related.
C                                    Note: this is a local index within
C                                    the x-type NBO atomic set, not the
C                                    overall atomic index of the entire
C                                    molecule.
C                    RING (A,B)   =  contains all atomic indices A of
C                                    those symmetry related sets that
C                                    have equal distance from center
C                                    B and are of maximum size.
C                    NRING (B)    =  contains the # of sets of maximum
C                                    size related to center B.
C                    RINGSZ (B)   =  contains the maximum size of the
C                                    above found sets related to
C                                    center B.
C                    SYMCEN       =  array that will be used to store
C                                    symmetry related atomic indices.
C                    LSYMACC      =  symmetry accuracy for orbital
C                                    localization contents
C                    QSYMACC      =  symmetry accuracy for orbital
C                                    interaction order values (sum of
C                                    squares of offdiagonal occupation
C                                    matrix elements for one orbital)
C                    NPAIR        =  local atomic index of the two-fold
C                                    weight degenerate atomic NBOs.
C                                    For each atom the local atomic
C                                    indices for the degenerate pair
C                                    are NPAIR and NPAIR+1.
C                    IDXDEG       =  array that will be used to store
C                                    the NBO indices that will be
C                                    rotated (for example, if sp2
C                                    hybrids will be constructed, this
C                                    array will indicate all the NBO
C                                    indices for the s and the two p
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
         LOGICAL     SAVEC,SAVEP
         LOGICAL     ZEROSYM

         INTEGER     AXIAL
         INTEGER     B1ST,BOND
         INTEGER     CENTER
         INTEGER     D,I,J,N,S
         INTEGER     LOBE
         INTEGER     MROT,NROT
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NBOSIZE,ROTSIZE
         INTEGER     NHO
         INTEGER     NLO__ANGULAR_MOMENTUM_NBO
         INTEGER     NPAIR
         INTEGER     NSYM
         INTEGER     NXA,NXBA
         INTEGER     OFFA,OFFX
         INTEGER     OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR
         INTEGER     POS
         INTEGER     SYMSIZE,SYMSETS
         INTEGER     TOTSYM

         INTEGER     ATXIDX  (1:NXA  )
         INTEGER     ATXOFF  (1:NXA  )
         INTEGER     BDNBAS  (1:NBAS )
         INTEGER     BDNCEN  (1:NBAS )
         INTEGER     IDXDEG  (1:NXBA )
         INTEGER     IDXSYM  (1:2    )
         INTEGER     NBOBD   (1:NBAS )
         INTEGER     NRING   (1:NATOM)
         INTEGER     RINGSZ  (1:NATOM)
         INTEGER     ROTIDX  (1:NBAS )
         INTEGER     SLTYPE  (1:NBAS )
         INTEGER     SYMCEN  (1:NATOM)
         INTEGER     SYMNBO  (1:NBAS )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS )
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS )
         INTEGER     RING    (1:NATOM  ,1:NATOM)

         DOUBLE PRECISION  HIJ
         DOUBLE PRECISION  LSYMACC,QSYMACC
         DOUBLE PRECISION  R11,R12,R21,R22
         DOUBLE PRECISION  SYMCRIT
         DOUBLE PRECISION  X,Y
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
C             ...set some data.
C
C
         SAVEC = .TRUE.
         SAVEP = .TRUE.
C
C
C             ...loop over both axial centers.
C
C
         DO AXIAL = 1,2

            CENTER = ATXIDX (IDXSYM (AXIAL))
            OFFA   = ATXOFF (IDXSYM (AXIAL))

            SYMSIZE = RINGSZ (CENTER)

            IF (SYMSIZE.LT.3) THEN
                WRITE (*,*) ' Symmetry conflict in axial case! '
                WRITE (*,*) ' CENTER,SYMSIZE = ',CENTER,SYMSIZE
                WRITE (*,*) ' nlo__rotate_axial_nbo '
                WRITE (1,*) ' Symmetry conflict in axial case! '
                WRITE (1,*) ' CENTER,SYMSIZE = ',CENTER,SYMSIZE
                WRITE (1,*) ' nlo__rotate_axial_nbo '
                STOP
            END IF

            SYMSETS = NRING (CENTER)
            NSYM = SYMSIZE * SYMSETS

            DO N = 1,NSYM
               SYMCEN (N) = RING (N,CENTER)
            END DO
C
C
C             ...rotate the degenerate pair first.
C
C
            D = OFFX + OFFA + NPAIR

            CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                 ( NBAS,NBAS,
     +                   NBAS,NBAS,
     +                   NBAS,2,
     +                   NBAS,2,
     +                   NBAS,
     +                   NBAS,NBAS,2,
     +                   0,0,
     +                   SAVEC,SAVEP,SAVEC,
     +                   C,P,C (1,D),
     +                   XVEC,
     +
     +                            XMAT )
     +
     +
            CALL  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                 ( NBAS,NATOM,2,
     +                   NBOSIZE,
     +                   OFFLP,OFFBD,OFFEP,
     +                   OFFAB,OFFRY,OFFRR,
     +                   CENTER,
     +                   SYMSIZE,SYMSETS,SYMCEN,
     +                   SYMNBO,
     +                   QSYMACC,LSYMACC,
     +                   BDNCEN,BDCEN,
     +                   NBOBD,
     +                   SYMCRIT,
     +                   XMAT,
     +                   LOCAL,
     +                   Q,
     +
     +                            ZEROSYM,
     +                            NROT,
     +                            ROTIDX )
     +
     +
            IF (ZEROSYM) THEN
                WRITE (*,*) ' NO SYMMETRY INTERACTION! '
                WRITE (*,*) ' Rotation failed! '
                RETURN
            END IF

            IF (NROT.EQ.0) THEN
                WRITE (*,*) ' Cannot establish NBO rotation set! '
                WRITE (*,*) ' CENTER = ',CENTER
                WRITE (*,*) ' nlo__rotate_axial_nbo '
                WRITE (1,*) ' Cannot establish NBO rotation set! '
                WRITE (1,*) ' CENTER = ',CENTER
                WRITE (1,*) ' nlo__rotate_axial_nbo '
                STOP
            END IF

            LOBE = ROTIDX (1)
            X    = XMAT (LOBE,1)
            Y    = XMAT (LOBE,2)

            MAXIMZE = .TRUE.

            CALL  NLO__FIND_2DEG_AXIAL_ROTMAT
     +
     +                 ( X,Y,
     +                   MAXIMZE,
     +
     +                           R11,R12,
     +                           R21,R22 )
     +
     +
            DO N = 1,NBAS
               X = R11 * C (N,D) + R21 * C (N,D+1)
               Y = R12 * C (N,D) + R22 * C (N,D+1)
               C (N,D  ) = X
               C (N,D+1) = Y
            END DO

            ROT (1,1) = R11
            ROT (2,1) = R21
            ROT (1,2) = R12
            ROT (2,2) = R22

            HYB (1,1) = ONE
            HYB (2,1) = ZERO
            HYB (1,2) = ZERO
            HYB (2,2) = ONE

            MROT = 2

            IDXDEG (1) = D
            IDXDEG (2) = D + 1
C
C
C             ...form the necessary hybrids, if possible and
C                necessary. Hybrids for 3-fold symmetry.
C
C
            IF (SYMSIZE.EQ.3) THEN

                IF (AXIAL.EQ.1) THEN

                    TOTSYM = NLO__ANGULAR_MOMENTUM_NBO
     +
     +                            ( NXBA,
     +                              SLTYPE (OFFX+OFFA+1),
     +                              NPAIR,NPAIR+1,
     +                              0,
     +                              1,
     +                              QSYMACC,LSYMACC,
     +                              W (OFFX+OFFA+1),
     +                              Q (OFFX+OFFA+1),
     +                              LOCAL (OFFX+OFFA+1) )
     +
     +
                    HYBRID = TOTSYM .NE. 0

                    IF (HYBRID .AND. ROTSIZE.LT.3) THEN
                        WRITE (*,*) ' Cannot form sp2 hybrids! '
                        WRITE (*,*) ' ROTSIZE = ',ROTSIZE
                        WRITE (*,*) ' nlo__rotate_axial_nbo '
                        WRITE (1,*) ' Cannot form sp2 hybrids! '
                        WRITE (1,*) ' ROTSIZE = ',ROTSIZE
                        WRITE (1,*) ' nlo__rotate_axial_nbo '
                        STOP
                    END IF
                END IF

                IF (HYBRID) THEN

                    S = OFFX + OFFA + TOTSYM

                    CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +                         ( ROTSIZE,ROTSIZE,
     +                           3,3,
     +
     +                                    ROT )
     +
     +

                    ROT (2,2) = R11
                    ROT (3,2) = R21
                    ROT (2,3) = R12
                    ROT (3,3) = R22

                    CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                         ( NBAS,NBAS,
     +                           NBAS,NBAS,
     +                           NBAS,1,
     +                           NBAS,1,
     +                           NBAS,
     +                           NBAS,NBAS,1,
     +                           0,0,
     +                           SAVEC,SAVEP,SAVEC,
     +                           C,P,C (1,S),
     +                           XVEC,
     +
     +                                    XMAT )
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

                    CALL  NLO__FORM_SP2_NBO
     +
     +                         ( NBAS,
     +
     +                                  HYB (1,1),HYB (1,2),HYB (1,3),
     +                                  HYB (2,1),HYB (2,2),HYB (2,3),
     +                                  HYB (3,1),HYB (3,2),HYB (3,3),
     +                                  C (1,S),
     +                                  C (1,D) )
     +
     +
                    MROT = 3

                    IDXDEG (1) = S
                    IDXDEG (2) = D
                    IDXDEG (3) = D + 1

                    SLTYPE (S  ) = -1
                    SLTYPE (D  ) = -1
                    SLTYPE (D+1) = -1

                ELSE

                    WRITE (*,*) ' No s function for sp2 hybrid! '
                    WRITE (*,*) ' Occuring at CENTER = ',CENTER
                    WRITE (*,*) ' Symmetry reduction C3v -> C2v! '
                    WRITE (1,*) ' No s function for sp2 hybrid! '
                    WRITE (1,*) ' Occuring at CENTER = ',CENTER
                    WRITE (1,*) ' Symmetry reduction C3v -> C2v! '

                END IF
C
C
C             ...No hybrids for 4-fold symmetry!.
C
C
            ELSE IF (SYMSIZE.GT.4) THEN

                WRITE (*,*) ' Cannot deal with SYMSIZE > 4! '
                WRITE (*,*) ' CENTER,SYMSIZE = ',CENTER,SYMSIZE
                WRITE (*,*) ' nlo__rotate_axial_nbo '
                WRITE (1,*) ' Cannot deal with SYMSIZE > 4! '
                WRITE (1,*) ' CENTER,SYMSIZE = ',CENTER,SYMSIZE
                WRITE (1,*) ' nlo__rotate_axial_nbo '
                STOP

            END IF
C
C
C             ...form the NBO coefficient matrix in NHO basis.
C
C
            DO J = 1,MROT
               POS = IDXDEG (J)
               DO N = 1,MROT
                  B (N,POS) = ZERO
               END DO
            END DO

            DO J = 1,MROT
               POS = IDXDEG (J)
               DO I = 1,MROT
                  HIJ = HYB (I,J)
                  DO N = 1,MROT
                     B (N,POS) = B (N,POS) + HIJ * ROT (N,I)
                  END DO
               END DO
            END DO
C
C
C             ...the matrices containing details about the new NBOs.
C
C
            POS = IDXDEG (1)
            B1ST = NBOBD (POS)
            BDNBAS (B1ST) = MROT

            DO I = 1,MROT
               POS = IDXDEG (I)
               BOND = NBOBD (POS)
               NHO = BDBAS (1,BOND)
               BDBAS (I,B1ST) = NHO
            END DO

            DO J = 2,MROT
               POS = IDXDEG (J)
               BOND = NBOBD (POS)
               BDNBAS (BOND) = MROT
               DO I = 1,MROT
                  BDBAS (I,BOND) = BDBAS (I,B1ST)
               END DO
            END DO
C
C
C             ...mark new NBOs as symmetry adapted.
C
C
            DO I = 1,MROT
               POS = IDXDEG (I)
               SYMNBO (POS) = 2
            END DO
C
C
C             ...next axial center.
C
C
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
