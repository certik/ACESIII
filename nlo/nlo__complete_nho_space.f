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
         SUBROUTINE  NLO__COMPLETE_NHO_SPACE
     +
     +                    ( NBAS,NHATOM,
     +                      NBOSIZE,
     +                      MXNHBA,MXNBA,
     +                      ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                      ATHIDX,
     +                      NHB,NHA,
     +                      NHYB,
     +                      WCRIT,WSTAR,
     +                      SAH,
     +                      P,PH,PHSUB,
     +                      W,C,
     +                      IVEC,
     +                      XVEC,XMAT,
     +
     +                              FAILED,
     +                              MORE,
     +                              NBOND,
     +                              BDNCEN,
     +                              BDCEN,
     +                              BDNBAS,
     +                              BDBAS,
     +                              BDOCC,
     +                              H )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__COMPLETE_NHO_SPACE
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine tries to complete the pre-NHO space on
C                each atomic center. The current bond forming pre-NHOs
C                on each atomic site have already been determined at
C                this stage, however NHATOM atoms still have an
C                incomplete bond forming pre-NHO set. Therefore, before
C                attempting to complete the atomic pre-NHO spaces, we
C                have to make sure that no more bonds can be formed
C                between all the NHATOM remaining centers. Two cases
C                can happen:
C
C                  1) If no more bonds can be formed, then the routine
C                     proceeds with completing the atomic pre-NHO
C                     spaces with possible Empty-pair bonds and Rydberg
C                     NHOs.
C
C                  2) If more bonds can still be formed, then we know
C                     that the latest bond and antibond delimiting
C                     weights have been not enough to find all bonds.
C                     In this case we ring the failure bell and exit.
C
C                For case 1) we proceed like follows:
C
C                  i) Extract the atomic submatrix PHSUB of the depleted
C                     PH occupation matrix and diagonalize it.
C
C                 ii) Starting with the ones corresponding to the
C                     largest eigenvalues, add as many eigenfunctions
C                     to the atomic pre-NHO set as needed to complete
C                     the bond forming pre-NHO set, i.e. until their
C                     number equals the number of Valence NAOs on that
C                     atom.
C
C                iii) WSW orthonormalize each atomic bond forming set of
C                     pre-NHOs to produce the bond forming NHOs.
C
C                 iv) Check the obtained Empty-pair NHOs as to what
C                     final weight they will have in the NBO density
C                     matrix. If any of the Empty-pair NHO weights
C                     is > WSTAR, we know that the present NHO
C                     formation pattern is inadequate. We issue the
C                     failure command and exit.
C
C                  v) Place +1 in the remaining atomic diagonal
C                     submatrix corresponding to the Rydberg NAOs and
C                     perform a partial Schmidt orthogonalization
C                     to the already orthonormal bond forming NHOs.
C
C                 vi) Normalize the orthogonal Rydberg NHOs obtained
C                     in v).
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NHATOM       =  current total # of atomic hybrid
C                                    centers left for checking.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    ATHCEN (I)   =  current hybrid atomic labels
C                                    (indices) for I-th atomic hybrid
C                                    center within the set of NHATOM
C                                    atomic hybrid centers still having
C                                    an incomplete bond forming pre-NHO
C                                    set.
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid NAOs
C                                    on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    ATHIDX (I)   =  will contain the reordered NHATOM
C                                    atomic indices correponding to
C                                    the unordered NHATOM atomic labels
C                                    present in array ATHCEN.
C                    NHB          =  total # of Hybrid NAOs.
C                    NHA          =  total # of hybrid atoms.
C                    NHYB (A)     =  current total # of hybrid NAO's
C                                    (pre-NHO's) on each atom A.
C                    WCRIT        =  abolute lowest weight above which a
C                                    bond formation is accepted.
C                    WSTAR        =  highest weight below which an
C                                    Empty-pair NHO formation is
C                                    accepted.
C                    SAH          =  will contain pre-NHO atomic overlap
C                                    matrices when checking for linear
C                                    dependencies on one atomic site.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    PH           =  current depleted NHB x NHB hybrid
C                                    occupation matrix.
C                    PHSUB        =  will contain the submatrices of
C                                    the occupation matrix PH.
C                    W            =  accumulated pre-NHO weight vector
C                                    for WSW procedure in the high
C                                    weight pre-NHO part.
C                    C            =  NAO coefficient matrix in AO basis
C                                    with columns in NHB/NCB/NRB order.
C                    IVEC         =  int scratch array of vector type
C                    XVEC         =  flp scratch array of vector type
C                    XMAT         =  flp scratch array of matrix type
C
C
C                  Output:
C
C                    FAILED       =  is set true, if the routine sees
C                                    an unrecoverable error in finishing
C                                    completion of the NHO sets.
C                    MORE         =  is true, if still at least one
C                                    more potential bond can be formed
C                                    between those atoms with incomplete
C                                    bond forming pre-NHO sets.
C                    NBOND        =  final total # of bonds.
C                    BDNCEN (J)   =  final # of atomic centers for
C                                    J-th bond.
C                    BDCEN (I,J)  =  final I-th atomic center index
C                                    for J-th bond.
C                    BDNBAS (J)   =  final # of basis functions (NHOs)
C                                    for J-th bond.
C                    BDBAS (I,J)  =  final I-th global basis (NHO)
C                                    index for J-th bond.
C                    BDOCC (I)    =  final # of occupied levels for
C                                    I-th bond.
C                    H (I,J)      =  MXNBA x NHB matrix containing the
C                                    final atomic NHOs corresponding
C                                    to the NHB space. I is the local
C                                    atomic index labeling the atomic
C                                    hybrid NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB space
C                                    NHOs, with all NHOs belonging to
C                                    a specific atomic center being
C                                    grouped together.
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

         LOGICAL     ABSOLUT
         LOGICAL     DEPEND
         LOGICAL     EXTRACT
         LOGICAL     FAILED
         LOGICAL     INCRESE
         LOGICAL     LOWDIN
         LOGICAL     LTRG,UTRG
         LOGICAL     MORE
         LOGICAL     NOOVLP
         LOGICAL     REVERS
         LOGICAL     SAVEH,SAVES,SAVEP

         INTEGER     ATOM
         INTEGER     I,J,K,L,M
         INTEGER     MXNHBA,MXNBA
         INTEGER     NBAS
         INTEGER     NBOND
         INTEGER     NBOSIZE
         INTEGER     NHATOM
         INTEGER     NHB,NHA
         INTEGER     NHBA,NEBA,NYBA
         INTEGER     NHOIDX
         INTEGER     NHYBA
         INTEGER     NVAL
         INTEGER     OFF

         INTEGER     ATHCEN  (1:NHA)
         INTEGER     ATHIDX  (1:NHATOM)
         INTEGER     ATNHB   (1:NHA)
         INTEGER     ATHOFF  (1:NHA)
         INTEGER     ATHVAL  (1:NHA)
         INTEGER     BDNBAS  (1:NHB)
         INTEGER     BDNCEN  (1:NHB)
         INTEGER     BDOCC   (1:NHB)
         INTEGER     IVEC    (1:NHATOM)
         INTEGER     NHYB    (1:NHA)

         INTEGER     BDBAS   (1:NBOSIZE,1:NHB)
         INTEGER     BDCEN   (1:NBOSIZE,1:NHB)

         DOUBLE PRECISION  E
         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  WCRIT,WSTAR,WEIGHT,WTHRESH
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  XVEC (1:2*NHB)
         DOUBLE PRECISION  W    (1:NHB)

         DOUBLE PRECISION  C     (1:NBAS  ,1:NBAS  )
         DOUBLE PRECISION  H     (1:MXNBA ,1:NHB   )
         DOUBLE PRECISION  P     (1:NBAS  ,1:NBAS  )
         DOUBLE PRECISION  PH    (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHSUB (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  SAH   (1:MXNHBA,1:MXNHBA)
         DOUBLE PRECISION  XMAT  (1:NBAS  ,1:MXNHBA)

         DATA  ONE      /1.D0/
         DATA  ZERO     /0.D0/
         DATA  REVERS   /.TRUE./
         DATA  WTHRESH  /1.D-4/
C
C
C------------------------------------------------------------------------
C
C
C             ...extract lower triangle of PHSUB matrix from PH matrix
C                corresponding to all hybrid atomic centers with
C                incomplete bond forming pre-NHO sets after reordering
C                of the atomic indices.
C
C
         ABSOLUT = .FALSE.
         INCRESE = .TRUE.

         CALL  NLO__SORT_INT_VECTOR_ELEMENTS
     +
     +              ( NHATOM,NHATOM,
     +                1,NHATOM,
     +                ABSOLUT,INCRESE,
     +                0,
     +                ATHCEN,
     +
     +                         IVEC )
     +
     +
         DO I = 1,NHATOM
            ATHIDX (I) = ATHCEN (IVEC (I))
         END DO

         LTRG = .TRUE.
         FAILED = .FALSE.
         EXTRACT = .TRUE.

         CALL  NLO__HANDLE_MATRIX_SECTIONS
     +
     +              ( NHB,NHB,
     +                NHB,NHB,
     +                NHA,NHATOM,
     +                NHATOM,
     +                ATHIDX,ATNHB,ATHOFF,
     +                EXTRACT,
     +                LTRG,
     +                PH,
     +
     +                        M,
     +                        PHSUB )
     +
     +
C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' PHSUB matrix for bond test',
C     +                  NHB,NHB,
C     +                  M,0,
C     +                  PHSUB )
C     +
C     +
C
C
C             ...diagonalize PHSUB matrix.
C
C
         CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +              ( NHB,NHB,NHB,
     +                M,
     +                REVERS,
     +
     +                        XVEC,
     +                        PHSUB )
     +
     +
C         CALL    MAT__PRINT_V_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' Eigenvalues of PHSUB matrix for bond test ',
C     +                  NHB,
C     +                  M,
C     +                  XVEC )
C     +
C     +
C
C
C             ...check eigenvalues for potential remaining bonds.
C                If at least one potential bond found, return.
C
C
         MORE = .FALSE.
         FAILED = .FALSE.

         DO I = 1,M
            IF (XVEC (I) .GT. WCRIT) THEN
                MORE = .TRUE.
                RETURN
            END IF
         END DO
C
C
C             ...start loop over all hybrid atomic sites, extract the
C                atomic submatrices PHSUB and analyze their eigenvalues.
C                The atomic overlap matrix over NAOs is equal to a unit
C                matrix, since the NAOs at this stage are orthonormal.
C                Add the Empty-pair and Rydberg bonds to the bond
C                characterization arrays.
C
C
         LTRG = .TRUE.
         UTRG  = .TRUE.
         SAVEH = .TRUE.
         SAVES = .TRUE.
         SAVEP = .TRUE.
         NOOVLP = .TRUE.
         EXTRACT = .TRUE.

         DO 1000 ATOM = 1,NHA

            OFF = ATHOFF (ATOM)
            NHBA = ATNHB (ATOM)
            NVAL = ATHVAL (ATOM)
            NHYBA = NHYB (ATOM)
            NYBA = NHBA - NVAL
            NEBA = NVAL - NHYBA
C
C
C             ...add the Empty-pair pre-NHOs (if necessary) to complete
C                the valence space.
C
C
            LOWDIN = .FALSE.

            IF (NEBA.GT.0) THEN

                CALL  NLO__HANDLE_MATRIX_SECTIONS
     +
     +                     ( NHB,NHB,
     +                       NHB,NHB,
     +                       NHA,1,
     +                       1,
     +                       ATOM,ATNHB,ATHOFF,
     +                       EXTRACT,
     +                       LTRG,
     +                       PH,
     +
     +                               M,
     +                               PHSUB )
     +
     +
                CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                     ( NHB,NHB,NHB,
     +                       M,
     +                       REVERS,
     +
     +                               XVEC,
     +                               PHSUB )
     +
     +
                K = NHYBA
                DO 100 J = 1,M
                   E = XVEC (J)
                   K = K + 1
                   NHOIDX = OFF + K
                   DO I = 1,M
                      H (I,NHOIDX) = PHSUB (I,J)
                   END DO

                   CALL  NLO__CHECK_LINEAR_DEPENDENCY
     +
     +                        ( MXNBA,K,
     +                          MXNHBA,MXNHBA,
     +                          M,K,
     +                          H (1,OFF+1),
     +                          SAH,
     +
     +                                 DEPEND )
     +
     +
                   IF (DEPEND) THEN
                       K = K - 1
                       GOTO 100
                   ELSE
                       W (NHOIDX) = MAX (E,WTHRESH)
                       NBOND = NBOND + 1
                       BDOCC (NBOND) = 0
                       BDNBAS (NBOND) = 1
                       BDNCEN (NBOND) = 1
                       BDCEN (1,NBOND) = ATOM
                       BDBAS (1,NBOND) = NHOIDX
                   END IF

                   IF (K.EQ.NVAL) GOTO 9000

  100           CONTINUE

                WRITE (*,*) ' Failure in adding pre-NHO Empty-pairs! '
                WRITE (*,*) ' Hybrid atom # = ',ATOM
                WRITE (*,*) ' nlo__complete_nho_space '
                WRITE (1,*) ' Failure in adding pre-NHO Empty-pairs! '
                WRITE (1,*) ' Hybrid atom # = ',ATOM
                WRITE (1,*) ' nlo__complete_nho_space '
                STOP

            END IF
C
C
C             ...perform the WSW orthonormalization on the valence
C                space.
C
C
 9000       CALL  NLO__WSW_ORTHONORMALIZE
     +
     +                 ( MXNBA,NVAL,
     +                   MXNHBA,MXNHBA,
     +                   NVAL,
     +                   NHBA,NVAL,
     +                   SAH,
     +                   W (OFF+1),
     +                   LOWDIN,
     +                   NOOVLP,
     +                   SAVES,
     +                   XVEC,XVEC (NVAL+1),
     +                   XMAT,
     +
     +                           FAILED,
     +                           H (1,OFF+1) )
     +
     +
            IF (FAILED) THEN
                WRITE (*,*) ' WSW orthonormalize failed for NVAL part! '
                WRITE (*,*) ' Hybrid atom # = ',ATOM
                WRITE (*,*) ' nlo__complete_nho_space '
                WRITE (1,*) ' WSW orthonormalize failed for NVAL part! '
                WRITE (1,*) ' Hybrid atom # = ',ATOM
                WRITE (1,*) ' nlo__complete_nho_space '
            END IF
C
C
C             ...check the obtained Empty-pair NHO weights.
C
C
            IF (NEBA.GT.0) THEN

                CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                     ( NBAS,NHBA,
     +                       NBAS,NBAS,
     +                       NBAS,NHBA,
     +                       NHBA,
     +                       NHBA,NBAS,
     +                       0,
     +                       SAVEP,LTRG,UTRG,
     +                       C (1,OFF+1),
     +                       P,
     +                       XVEC,
     +
     +                              XMAT )
     +
     +
                DO J = 1,NEBA
                   NHOIDX = OFF + NHYBA + J
                   WEIGHT = ZERO
                   DO L = 1,NHBA
                      SUM = ZERO
                      DO K = 1,NHBA
                         SUM = SUM + H (K,NHOIDX) * XMAT (K,L)
                      END DO
                      WEIGHT = WEIGHT + SUM * H (L,NHOIDX)
                   END DO

                   IF (WEIGHT.GT.WSTAR) THEN
                       WRITE (*,*) ' Empty pair weight > WSTAR! '
                       WRITE (*,*) ' Hybrid atom # = ',ATOM
                       WRITE (*,*) ' nlo__complete_nho_space '
                       WRITE (1,*) ' Empty pair weight > WSTAR! '
                       WRITE (1,*) ' Hybrid atom # = ',ATOM
                       WRITE (1,*) ' nlo__complete_nho_space '
                       FAILED = .TRUE.
                       RETURN
                   END IF

                END DO

            END IF
C
C
C             ...add the Rydberg pre-NHOs (if necessary) to complete
C                the entire NHO space on the present hybrid atom and
C                perform a partial Schmidt orthogonalization of the
C                Rydberg pre-NHOs to the already orthogonal Valence
C                NHOs plus a subsequent Lowdin orthonormalization of
C                the Rydberg pre-NHOs. The set of orthonormal Rydberg
C                NHOs thus obtained has to be further linearly
C                combined to reflect the right molecular symmetry.
C                To this end the occupation matrix is set up in the
C                Rydberg NHO space and diagonalized. The eigenvectors
C                define the final Rydberg NHOs.
C
C
            LOWDIN = .TRUE.

            IF (NYBA.GT.0) THEN

                DO J = 1,NYBA
                   NHOIDX = OFF + NVAL + J
                   H (NVAL+J,NHOIDX) = ONE
                   NBOND = NBOND + 1
                   BDOCC (NBOND) = 0
                   BDNBAS (NBOND) = 1
                   BDNCEN (NBOND) = 1
                   BDCEN (1,NBOND) = ATOM
                   BDBAS (1,NBOND) = NHOIDX
                END DO

                CALL  MAT__C_EQ_UNIT_FLOAT
     +
     +                     ( MXNHBA,MXNHBA,
     +                       MXNHBA,MXNHBA,
     +
     +                                SAH )
     +
     +
                CALL  NLO__PARTIAL_SCHMIDT
     +
     +                     ( MXNBA,NVAL,
     +                       MXNBA,NYBA,
     +                       MXNHBA,MXNHBA,
     +                       NHBA,NVAL,NYBA,
     +                       MAX (NVAL,NYBA),MIN (NVAL,NYBA),
     +                       SAH,
     +                       SAVES,SAVEH,SAVEH,
     +                       H (1,OFF+1),
     +                       XVEC,
     +                       XMAT,
     +
     +                                H (1,OFF+NVAL+1) )
     +
     +
                CALL  NLO__WSW_ORTHONORMALIZE
     +
     +                     ( MXNBA,NYBA,
     +                       MXNHBA,MXNHBA,
     +                       NYBA,
     +                       NHBA,NYBA,
     +                       SAH,
     +                       W (OFF+NVAL+1),
     +                       LOWDIN,
     +                       NOOVLP,
     +                       SAVES,
     +                       XVEC,XVEC (NYBA+1),
     +                       XMAT,
     +
     +                                FAILED,
     +                                H (1,OFF+NVAL+1) )
     +
     +
                IF (FAILED) THEN
                    WRITE (*,*) ' WSW process failed at Rydberg part! '
                    WRITE (*,*) ' Hybrid atom # = ',ATOM
                    WRITE (*,*) ' nlo__complete_nho_space '
                    WRITE (1,*) ' WSW process failed at Rydberg part! '
                    WRITE (1,*) ' Hybrid atom # = ',ATOM
                    WRITE (1,*) ' nlo__complete_nho_space '
                    RETURN
                END IF

                CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                     ( NBAS,NHBA,
     +                       NBAS,NBAS,
     +                       NBAS,NHBA,
     +                       NHBA,
     +                       NHBA,NBAS,
     +                       0,
     +                       SAVEP,LTRG,UTRG,
     +                       C (1,OFF+1),
     +                       P,
     +                       XVEC,
     +
     +                              XMAT )
     +
     +
                CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                     ( MXNBA,NYBA,
     +                       NBAS,NHBA,
     +                       MXNHBA,NYBA,
     +                       NYBA,
     +                       NYBA,NHBA,
     +                       0,
     +                       SAVEP,LTRG,UTRG,
     +                       H (1,OFF+NVAL+1),
     +                       XMAT,
     +                       XVEC,
     +
     +                              SAH )
     +
     +
                CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                     ( MXNHBA,NYBA,NYBA,
     +                       NYBA,
     +                       REVERS,
     +
     +                               XVEC,
     +                               SAH )
     +
     +
                CALL  MAT__C_EQ_A_FLOAT
     +
     +                     ( MXNBA,NYBA,
     +                       NBAS,NYBA,
     +                       NHBA,NYBA,
     +                       H (1,OFF+NVAL+1),
     +
     +                               XMAT )
     +
     +
                CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                     ( NBAS,NYBA,
     +                       MXNHBA,NYBA,
     +                       MXNBA,NYBA,
     +                       NHBA,NYBA,NYBA,
     +                       XMAT,SAH,
     +
     +                               H (1,OFF+NVAL+1) )
     +
     +
            END IF
C
C
C             ...next hybrid atomic site.
C
C
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
