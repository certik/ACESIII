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
         SUBROUTINE  NLO__GENER_NAO_ORBITALS
     +
     +                    ( NBAS,NATOM,
     +                      MXSHELL,MXLSIZE,MXNAL,MXNBA,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      P,PA,S,
     +                      SAL,CAL,WAL,
     +                      TSYM,
     +                      MJUMP,
     +                      NMB,NRB,
     +                      LSIZE,
     +                      BASBEG,BASEND,
     +                      PSYMACC,
     +                      SYMMAP,
     +                      ATDONE,IDXSYM,
     +                      COLMAP,
     +                      NRYDAL,
     +                      IVEC,XVEC,XMAT,
     +                      WPRE,
     +
     +                              BDSIZE,
     +                              BDATOM,
     +                              ANGSYM,
     +                              W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NAO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the final set of NAO's in
C                atomic order from the input set of pre-NAO's in atomic
C                order. The steps performed are:
C
C
C                   1) Reorder the weights and pre-NAO coefficients
C                      from atomic to NMB/NRB order.
C
C
C                        ----- Orthonormalization of pre-NAO set ------
C
C
C                   2) Peform a WSW (weighted) orthonormalization on
C                      the NMB pre-NAO's to obtain the orthonormal
C                      NMB pre-NAO's.
C
C                           If NRB = 0, jump directly to step 10)
C
C                   3) Schmidt orthogonalize the NRB pre-NAO's to the
C                      orthogonal NMB pre-NAO's.
C
C                   4) Restore the pre-NAO character of the NRB space
C                      by diagonalizing the systems of equations
C                      P(AL)C(AL) = S(AL)C(AL)W(AL) for all Rydberg
C                      (AL) spaces.
C
C                   5) Divide the NRB pre-NAO's into two sets NRBINT
C                      and NRBEXT, depending on their weights.
C
C                   6) Reorder the NRB weigths and pre-NAO coefficients
C                      from NRB to NRBINT/NRBEXT order.
C
C                   7) Perform a scaled WSW (weighted) orthonormalization
C                      on the NRBINT pre-NAO's to obtain new orthonormal
C                      NRBINT pre-NAO's.
C
C                   8) Schmidt orthogonalize the NRBEXT pre-NAO's to the
C                      orthonormal NRBINT pre-NAO's.
C
C                   9) Perform a Loewdin orthonormalization on the
C                      NRBEXT pre-NAO's to obtain new orthonormal
C                      NRBEXT NAO's.
C
C                  10) Reorder the orthonormal pre-NAO coefficients from
C                      NMB/NRBINT/NRBEXT to atomic order.
C
C
C                        ----- Formation of final NAO set ------
C
C
C                  11) Form the NAO's by diagonalizing the systems of
C                      equations P(AL)C(AL) = S(AL)C(AL)W(AL) for all
C                      (AL) spaces.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MXLSIZE      =  largest l-size value, i.e. maximum
C                                    m-degeneracy.
C                    MXNAL        =  maximum size of atomic l-shell
C                                    space. The atomic l-shell space
C                                    is the total # of contractions for
C                                    an atomic l-shell.
C                    MXNBA        =  maximum # of AOs per atom
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A (in increasing order!).
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    PA           =  atomic occupation matrix.
C                    S            =  full NBAS x NBAS overlap matrix.
C                    SAL,CAL,WAL  =  submatrices S(AL) and C(AL) for
C                                    overlap and density/pre-NAO coeffs
C                                    and subvector W(AL) for weights.
C                    TSYM         =  MXLSIZE x MXLSIZE matrix that will
C                                    be used to hold the l-sitesymmetry
C                                    transformation matrix.
C                    MJUMP        =  is .true., if the m values in the
C                                    m-space are ordered such that the
C                                    same m values are separated. This
C                                    keyword is necessary because some
C                                    AO basis functions are m-ordered
C                                    differently within each l-shell.
C                                    It invokes different types of
C                                    m-averaging algorithms.
C                    NMB,NRB      =  # of minimal and Rydberg type
C                                    pre-NAO's.
C                    LSIZE (I)    =  I-th l-shell size
C                    BASBEG (A)   =  first basis index number for atom A
C                    BASEND (A)   =  last basis index number for atom A
C                    PSYMACC      =  symmetry accuracy for occupation
C                                    matrix values (including weights)
C                    SYMMAP (A,B) =  has been set equal to 1 if atoms
C                                    A and B were found to be symmetry
C                                    related at the present stage of
C                                    calculation. A value of 0 indicates
C                                    no symmetry relation up to now.
C                    ATDONE (A)   =  will be used as an indicator if
C                                    atom A has been checked for
C                                    symmetry or not.
C                    IDXSYM       =  array that will be used to store
C                                    symmetry related atomic indices.
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of I-th pre-NAO,
C                                    based on atomic order, in the
C                                    NMB/NRB order.
C                    NRYDAL (I,A) =  # of Rydberg pre-NAO's found for
C                                    the I-th atomic l-shell for atom A.
C                    IVEC,XVEC    =  int/flp scratch array of vector
C                                    type.
C                    XMAT         =  flp scratch array of matrix type.
C                    WPRE         =  pre-NAO weight vector in atomic
C                                    order.
C                    C            =  pre-NAO coefficient matrix in AO
C                                    basis with columns in atomic order.
C
C
C                  Output:
C
C                    BDSIZE (I)   =  NAO bond sizes, i.e. the # of atoms
C                                    which form the I-th NAO. Of course,
C                                    this is trivially = 1 for all NAO. 
C                    BDATOM (I)   =  Atomic index map for the NAOs.
C                                    Atomic index forming the I-th
C                                    NAO.
C                    ANGSYM (I)   =  angular symmetry of I-th NAO.
C                                    This is equal to the lowest l-shell
C                                    value that can mix with the I-th
C                                    NAO and is dependent on the site
C                                    symmetry of the atom to which the
C                                    I-th NAO belongs.
C                    W            =  NAO weight vector in atomic order.
C                    C            =  NAO coefficient matrix in AO basis
C                                    with columns in atomic order.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ADD
         LOGICAL     FAILED
         LOGICAL     LOWDIN
         LOGICAL     LTRG,UTRG
         LOGICAL     MJUMP
         LOGICAL     NEW
         LOGICAL     NOOVLP
         LOGICAL     ORDER
         LOGICAL     REVERS
         LOGICAL     SAVEC,SAVEP,SAVES

         INTEGER     A,B
         INTEGER     ATOM
         INTEGER     BASB,BASL,BASNR
         INTEGER     I,J,K,N
         INTEGER     INDEX
         INTEGER     LDIM,LTOT,LTYPE
         INTEGER     MXSHELL,MXLSIZE,MXNAL,MXNBA
         INTEGER     NBAS,NATOM
         INTEGER     NAL
         INTEGER     NLTYPE
         INTEGER     NMB,NRB,NRBINT,NRBEXT
         INTEGER     NSYM
         INTEGER     NXBA

         INTEGER     ANGSYM  (1:NBAS     )
         INTEGER     ATDONE  (1:NATOM    )
         INTEGER     BASBEG  (1:NATOM    )
         INTEGER     BASEND  (1:NATOM    )
         INTEGER     BDATOM  (1:NBAS     )
         INTEGER     BDSIZE  (1:NBAS     )
         INTEGER     COLMAP  (1:NBAS     )
         INTEGER     IDXSYM  (1:NATOM    )
         INTEGER     IVEC    (1:NBAS+NBAS)
         INTEGER     LSIZE   (0:MXSHELL  )
         INTEGER     NSHELLS (1:NATOM    )

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     NRYDAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM    ,1:NATOM)

         DOUBLE PRECISION  PSYMACC
         DOUBLE PRECISION  WEIGHT,WTHRESH,WSCALE
         DOUBLE PRECISION  X
         DOUBLE PRECISION  ZERO,ONE,WEQZERO

         DOUBLE PRECISION  W    (1:NBAS     )
         DOUBLE PRECISION  WAL  (1:MXNAL    )
         DOUBLE PRECISION  WPRE (1:NBAS     )
         DOUBLE PRECISION  XVEC (1:NBAS+NBAS)

         DOUBLE PRECISION  C    (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  S    (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  P    (1:NBAS   ,1:NBAS   )
         DOUBLE PRECISION  PA   (1:MXNBA  ,1:MXNBA  )
         DOUBLE PRECISION  SAL  (1:MXNAL  ,1:MXNAL  )
         DOUBLE PRECISION  CAL  (1:MXNAL  ,1:MXNAL  )
         DOUBLE PRECISION  TSYM (1:MXLSIZE,1:MXLSIZE)
         DOUBLE PRECISION  XMAT (1:NBAS   ,1:NBAS   )

         DATA  ONE      /1.D0/
         DATA  ZERO     /0.D0/
         DATA  REVERS   /.TRUE./
         DATA  WTHRESH  /1.D-4/
         DATA  WEQZERO  /1.D-12/
C
C
C------------------------------------------------------------------------
C
C
C             ...save original pre-NAO weight vector in atomic order,
C                and reorder the pre-NAO coefficient matrix and the
C                weight vector from atomic to NMB/NRB order.
C
C
C         CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
C     +
C     +                ( 6,
C     +                  ' pre-NAO weight vector ',
C     +                  NBAS,
C     +                  NBAS,
C     +                  WPRE )
C     +
C     +
         CALL  MAT__W_EQ_U_FLOAT
     +
     +              ( NBAS,
     +                NBAS,
     +                NBAS,
     +                WPRE,
     +
     +                        W )
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
C             ...perform a WSW orthonormalization on the NMB part.
C
C
         SAVES = .TRUE.
         LOWDIN = .FALSE.
         NOOVLP = .FALSE.

         CALL    NLO__WSW_ORTHONORMALIZE
     +
     +                ( NBAS,NMB,
     +                  NBAS,NBAS,
     +                  NMB,
     +                  NBAS,NMB,
     +                  S,
     +                  W,
     +                  LOWDIN,
     +                  NOOVLP,
     +                  SAVES,
     +                  XVEC,XVEC (NMB+1),
     +                  XMAT,
     +
     +                          FAILED,
     +                          C )
     +
     +
         IF (FAILED) THEN
             WRITE (*,*) ' WSW orthonormalize failed for NMB part! '
             WRITE (*,*) ' nlo__gener_nao_orbitals '
             WRITE (1,*) ' WSW orthonormalize failed for NMB part! '
             WRITE (1,*) ' nlo__gener_nao_orbitals '
             STOP
         END IF
C
C
C             ...enter the Rydberg space processing (only, if Rydberg
C                space exists). Schmidt orthogonalize the NRB set of
C                pre-NAO's to the already orthonormal NMB set of
C                pre-NAO's.
C
C
         IF (NRB.NE.0) THEN

             SAVEC = .TRUE.
             SAVES = .TRUE.

             CALL    NLO__PARTIAL_SCHMIDT
     +
     +                    ( NBAS,NMB,
     +                      NBAS,NRB,
     +                      NBAS,NBAS,
     +                      NBAS,NMB,NRB,
     +                      MAX (NMB,NRB),MIN (NMB,NRB),
     +                      S,
     +                      SAVEC,SAVES,SAVEC,
     +                      C (1,1),
     +                      XVEC,
     +                      XMAT,
     +
     +                               C (1,NMB+1) )
     +
     +
C
C
C             ...restore the natural character of the Rydberg pre-NAO
C                basis, which has been modified during the above
C                Schmidt orthogonalization of the NRB pre-NAO set to
C                the orthogonal NMB pre-NAO set.
C
C                The structure of the algorithm is the same as for
C                the determination of the pre-NAO's, only this time
C                we look at the Rydberg space only. Symmetry related
C                atoms are treated on an equal footing.
C
C
             LTRG  = .TRUE.
             UTRG  = .FALSE.
             SAVEP = .TRUE.
             SAVES = .TRUE.

             BASB = NMB
             DO B = 1,NATOM
                IVEC (B) = BASB
                NLTYPE = NSHELLS (B)
                DO N = 1,NLTYPE
                   BASB = BASB + NRYDAL (N,B)
                END DO
             END DO

             DO A = 1,NATOM
                ATDONE (A) = 0
             END DO

             DO A = 1,NATOM
                NEW = ATDONE (A) .EQ. 0
                IF (NEW) THEN

                    NSYM = 0
                    DO B = 1,NATOM
                       IF (SYMMAP (B,A) .EQ. 1) THEN
                           NSYM = NSYM + 1
                           IDXSYM (NSYM) = B
                           ATDONE (B) = 1
                       END IF
                    END DO

                    BASL = 0
                    ATOM = IDXSYM (1)
                    NLTYPE = NSHELLS (ATOM)

                    DO N = 1,NLTYPE
                       LTOT = NRYDAL (N,ATOM)

                       IF (LTOT.NE.0) THEN
                           LTYPE = SHELLS (N,ATOM)
                           LDIM = LSIZE (LTYPE)
                           IF ( MOD (LTOT,LDIM).NE.0 ) THEN
                                WRITE (*,*) ' Bad NRB space size! '
                                WRITE (*,*) ' ATOM,SHELL = ',ATOM,LTYPE
                                WRITE (*,*) ' nlo__gener_nao_orbitals '
                                WRITE (1,*) ' Bad NRB space size! '
                                WRITE (1,*) ' ATOM,SHELL = ',ATOM,LTYPE
                                WRITE (1,*) ' nlo__gener_nao_orbitals '
                                STOP
                           END IF
                           NAL = LTOT / LDIM

                           DO I = 1,NSYM
                              B = IDXSYM (I)
                              BASB = IVEC (B)
                              BASNR = BASB + BASL
                              ADD = I.GT.1

                              CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                                   ( NBAS,LTOT,
     +                                     NBAS,NBAS,
     +                                     NBAS,LTOT,
     +                                     LTOT,
     +                                     LTOT,NBAS,
     +                                     0,
     +                                     SAVEP,LTRG,UTRG,
     +                                     C (1,BASNR+1),
     +                                     P,
     +                                     XVEC,
     +
     +                                             XMAT )
     +
     +
                              CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                                   ( NBAS,NBAS,
     +                                     MXNAL,MXNAL,
     +                                     LTOT,LDIM,
     +                                     0,
     +                                     MJUMP,LTRG,ADD,
     +                                     XMAT,
     +
     +                                             CAL )
     +
     +
                              CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                                   ( NBAS,LTOT,
     +                                     NBAS,NBAS,
     +                                     NBAS,LTOT,
     +                                     LTOT,
     +                                     LTOT,NBAS,
     +                                     0,
     +                                     SAVES,LTRG,UTRG,
     +                                     C (1,BASNR+1),
     +                                     S,
     +                                     XVEC,
     +
     +                                             XMAT )
     +
     +
                              CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                                   ( NBAS,NBAS,
     +                                     MXNAL,MXNAL,
     +                                     LTOT,LDIM,
     +                                     0,
     +                                     MJUMP,LTRG,ADD,
     +                                     XMAT,
     +
     +                                             SAL )
     +
     +
                           END DO

                           IF (NSYM.GT.1) THEN
                               X = ONE / DFLOAT (NSYM)
                               DO J = 1,NAL
                               DO I = J,NAL
                                  CAL (I,J) = X * CAL (I,J)
                                  SAL (I,J) = X * SAL (I,J)
                               END DO
                               END DO
                           END IF

                           CALL  MAT__GEN_EIGSYS_REAL_SYMMETRIC
     +
     +                                ( MXNAL,MXNAL,
     +                                  MXNAL,MXNAL,
     +                                  MXNAL,
     +                                  NAL,
     +                                  REVERS,
     +
     +                                          WAL,
     +                                          SAL,
     +                                          CAL )
     +
     +
                           DO I = 1,NSYM
                              B = IDXSYM (I)
                              BASB = IVEC (B)
                              BASNR = BASB + BASL

                              CALL  MAT__C_EQ_A_FLOAT
     +
     +                                   ( NBAS,LTOT,
     +                                     NBAS,LTOT,
     +                                     NBAS,LTOT,
     +                                     C (1,BASNR+1),
     +
     +                                             XMAT )
     +
     +
                              CALL  NLO__FORM_M_EXPANDED_COEFFS
     +
     +                                   ( NBAS,LTOT,
     +                                     MXNAL,MXNAL,
     +                                     NBAS,LTOT,
     +                                     NBAS,LTOT,NAL,
     +                                     MJUMP,
     +                                     XMAT,CAL,
     +
     +                                             C (1,BASNR+1) )
     +
     +
                              CALL  NLO__FORM_M_EXPANDED_WEIGHTS
     +
     +                                   ( MXNAL,
     +                                     LTOT,
     +                                     LTOT,NAL,
     +                                     MJUMP,
     +                                     WAL,
     +
     +                                             W (BASNR+1) )
     +
     +
                           END DO

                           BASL = BASL + LTOT

                       END IF
                    END DO
                END IF
             END DO
C
C
C             ...split the NRB Rydberg pre-NAO's into two sets
C                according to their weights. The first set will
C                contain those Rydberg pre-NAO's whose rescaled
C                weights (defined such that the largest weight
C                within the entire NRB set is equal to 1) are
C                larger than a threshold value WTHRESH. The second
C                set contains all the rest. The first set containing
C                NRBINT elements, will be subjected to a weighted
C                orthonormalization procedure with rescaled weights.
C                The rest, containing NRBEXT elements, will be
C                subjected to a normal Loewdin orthonormalization,
C                which can be performed with the same routine used
C                for weighted orthonormalizations but with all weights
C                equal to a constant value.
C
C                The first steps that follow are the rescaling of the
C                NRB weights, the decomposition of the NRB space into
C                NRBINT and NRBEXT, the determination of the Rydberg
C                column mapping (this is a local mapping!), the
C                local reordering of the Rydberg part of the pre-NAO
C                coefficient matrix and the update of the global
C                column permutation map.
C
C
             WSCALE = ZERO
             DO I = 1,NRB
                WSCALE = MAX (WSCALE,W (NMB+I))
             END DO
             WSCALE = ONE / WSCALE

             NRBINT = 0
             NRBEXT = 0

             DO I = 1,NRB
                WEIGHT = WSCALE * W (NMB+I)
                IF (WEIGHT.GT.WTHRESH) THEN
                    NRBINT = NRBINT + 1
                    IVEC (I) = NRBINT
                    XVEC (I) = WEIGHT
                ELSE
                    NRBEXT = NRBEXT + 1
                    IVEC (I) = NRBEXT + NBAS
                END IF
             END DO
         
             DO I = 1,NRB
                INDEX = IVEC (I)
                IF (INDEX.GT.NBAS) THEN
                    INDEX = INDEX - NBAS + NRBINT
                    IVEC (I) = INDEX
                ELSE
                    W (NMB+INDEX) = XVEC (I)
                END IF
             END DO

             CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +                  ( NBAS,NRB,
     +                    NRB,
     +                    NRB,
     +                    NBAS,
     +                    NBAS,NRB,
     +                    IVEC,
     +                    IVEC (NRB+1),
     +                    XVEC,
     +
     +                            C (1,NMB+1))
     +
     +
             DO I = 1,NBAS
                INDEX = COLMAP (I)
                IF (INDEX.GT.NMB) THEN
                    COLMAP (I) = IVEC (INDEX-NMB) + NMB
                END IF
             END DO
C
C
C             ...perform a WSW orthonormalization on the NRBINT part.
C
C
             SAVES = .TRUE.
             LOWDIN = .FALSE.
             NOOVLP = .FALSE.

             CALL  NLO__WSW_ORTHONORMALIZE
     +
     +                  ( NBAS,NRBINT,
     +                    NBAS,NBAS,
     +                    NRBINT,
     +                    NBAS,NRBINT,
     +                    S,
     +                    W (NMB+1),
     +                    LOWDIN,
     +                    NOOVLP,
     +                    SAVES,
     +                    XVEC,XVEC (NRBINT+1),
     +                    XMAT,
     +
     +                            FAILED,
     +                            C (1,NMB+1) )
     +
     +
             IF (FAILED) THEN
                 WRITE (*,*) ' WSW orthonorm failed for NRBINT part! '
                 WRITE (*,*) ' nlo__gener_nao_orbitals '
                 WRITE (1,*) ' WSW orthonorm failed for NRBINT part! '
                 WRITE (1,*) ' nlo__gener_nao_orbitals '
                 STOP
             END IF
C
C
C             ...Schmidt orthogonalize the NRBEXT set of pre-NAO's to
C                the already orthonormal NRBINT set of pre-NAO's.
C
C
             SAVEC = .TRUE.
             SAVES = .TRUE.

             CALL  NLO__PARTIAL_SCHMIDT
     +
     +                  ( NBAS,NRBINT,
     +                    NBAS,NRBEXT,
     +                    NBAS,NBAS,
     +                    NBAS,NRBINT,NRBEXT,
     +                    MAX (NRBINT,NRBEXT),MIN (NRBINT,NRBEXT),
     +                    S,
     +                    SAVEC,SAVES,SAVEC,
     +                    C (1,NMB+1),
     +                    XVEC,
     +                    XMAT,
     +
     +                             C (1,NMB+NRBINT+1) )
     +
     +
C
C
C             ...perform a WSW orthonormalization with equal weights
C                (i.e. a Loewdin orthonormalization) on the NRBEXT part.
C
C
             SAVES = .TRUE.
             LOWDIN = .TRUE.
             NOOVLP = .FALSE.

             CALL  NLO__WSW_ORTHONORMALIZE
     +
     +                  ( NBAS,NRBEXT,
     +                    NBAS,NBAS,
     +                    NRBEXT,
     +                    NBAS,NRBEXT,
     +                    S,
     +                    W (NMB+NRBINT+1),
     +                    LOWDIN,
     +                    NOOVLP,
     +                    SAVES,
     +                    XVEC,XVEC (NRBEXT+1),
     +                    XMAT,
     +
     +                            FAILED,
     +                            C (1,NMB+NRBINT+1) )
     +
     +
             IF (FAILED) THEN
                 WRITE (*,*) ' WSW orthonorm failed for NRBEXT part! '
                 WRITE (*,*) ' nlo__gener_nao_orbitals '
                 WRITE (1,*) ' WSW orthonorm failed for NRBEXT part! '
                 WRITE (1,*) ' nlo__gener_nao_orbitals '
                 STOP
             END IF
C
C
C             ...end Rydberg space manipulations.
C
C
         END IF
C
C
C             ...find invers of map: atomic -> NMB/NRBINT/NRBEXT
C                and restore the now orthonormal pre-NAO coefficient
C                matrix to atomic order.
C
C
         DO I = 1,NBAS
            IVEC (COLMAP (I)) = I
         END DO

         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,NBAS,
     +                IVEC,
     +                IVEC (NBAS+1),
     +                XVEC,
     +
     +                        C )
     +
     +
C
C
C             ...form the final NAO weights and coefficients in
C                atomic order. Note, that the present set of pre-NAO's
C                is now orthonormal, hence no atomic overlap submatrix
C                evaluation is necessary. Check the resulting NAO
C                weights for zeros to within numerical accuracy and
C                set them exactly equal to zero. This step is necessary
C                to avoid a later preliminary (unnecessary) stop
C                of the program when analyzing the NAO weights.
C                Set the NAO bond size array trivially equal to 1 and
C                fill in the NAO bond atom index array with the
C                appropriate atomic indices.
C
C                The structure of the code is essentially the same as
C                the one used to generate the initital pre-NAOs.
C                Symmetry equivalent atoms are treated together to
C                ensure NAO coefficient equalities.
C
C                Also, within each l-space, the NAOs will be ordered
C                such that the m-space is in separate order. For
C                example the p-shell NAOs will be ordered in the
C                form px,py,pz,px,py,pz,...This is achieved by
C                determining the column reordering map within each
C                l-space if MJUMP is .false.
C
C
         LTRG  = .TRUE.
         UTRG  = .FALSE.
         SAVEP = .TRUE.
         ORDER = .NOT.MJUMP

         DO A = 1,NATOM
            ATDONE (A) = 0
         END DO

         DO A = 1,NATOM
            NEW = ATDONE (A) .EQ. 0
            IF (NEW) THEN

                NSYM = 0
                DO B = 1,NATOM
                   IF (SYMMAP (B,A) .EQ. 1) THEN
                       NSYM = NSYM + 1
                       IDXSYM (NSYM) = B
                       ATDONE (B) = 1
                   END IF
                END DO

                BASL = 0
                ATOM = IDXSYM (1)
                NLTYPE = NSHELLS (ATOM)

                DO N = 1,NLTYPE
                   NAL = NBASAL (N,ATOM)
                   LTYPE = SHELLS (N,ATOM)
                   LDIM = LSIZE (LTYPE)
                   LTOT = NAL * LDIM

                   IF (ORDER) THEN
                       K = 0
                       DO I = 1,NAL
                       DO J = 1,LDIM
                          K = K + 1
                          COLMAP (K) = (J-1) * LDIM + I
                       END DO
                       END DO
                   END IF

                   DO I = 1,NSYM
                      B = IDXSYM (I)
                      BASB = BASBEG (B) - 1
                      BASNR = BASB + BASL
                      ADD = I.GT.1

                      CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                           ( NBAS,LTOT,
     +                             NBAS,NBAS,
     +                             NBAS,LTOT,
     +                             LTOT,
     +                             LTOT,NBAS,
     +                             0,
     +                             SAVEP,LTRG,UTRG,
     +                             C (1,BASNR+1),
     +                             P,
     +                             XVEC,
     +
     +                                     XMAT )
     +
     +
                      CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                           ( NBAS,NBAS,
     +                             MXNAL,MXNAL,
     +                             LTOT,LDIM,
     +                             0,
     +                             MJUMP,LTRG,ADD,
     +                             XMAT,
     +
     +                                     CAL )
     +
     +
                   END DO

                   IF (NSYM.GT.1) THEN
                       X = ONE / DFLOAT (NSYM)
                       DO J = 1,NAL
                       DO I = J,NAL
                          CAL (I,J) = X * CAL (I,J)
                       END DO
                       END DO
                   END IF

              WRITE (*,*) ' LTYPE,NAL = ',LTYPE,NAL

         CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +              ( 6,
     +                ' CAL matrix ',
     +                MXNAL,MXNAL,
     +                NAL,0,
     +                CAL )
     +
     +
                   CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                        ( MXNAL,MXNAL,MXNAL,
     +                          NAL,
     +                          REVERS,
     +
     +                                  WAL,
     +                                  CAL )
     +
     +
         CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +              ( 6,
     +                ' CAL eigenvalues ',
     +                MXNAL,1,
     +                NAL,1,
     +                WAL )
     +
     +
         CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +              ( 6,
     +                ' CAL eigenvectors ',
     +                MXNAL,MXNAL,
     +                NAL,NAL,
     +                CAL )
     +
     +
                   DO I = 1,NSYM
                      B = IDXSYM (I)
                      BASB = BASBEG (B) - 1
                      BASNR = BASB + BASL

                      CALL  MAT__C_EQ_A_FLOAT
     +
     +                           ( NBAS,LTOT,
     +                             NBAS,LTOT,
     +                             NBAS,LTOT,
     +                             C (1,BASNR+1),
     +
     +                                     XMAT )
     +
     +
                      CALL  NLO__FORM_M_EXPANDED_COEFFS
     +
     +                           ( NBAS,LTOT,
     +                             MXNAL,MXNAL,
     +                             NBAS,LTOT,
     +                             NBAS,LTOT,NAL,
     +                             MJUMP,
     +                             XMAT,CAL,
     +
     +                                     C (1,BASNR+1) )
     +
     +
                      IF (ORDER) THEN
                          CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +                               ( NBAS,LTOT,
     +                                 LTOT,
     +                                 LTOT,
     +                                 NBAS,
     +                                 NBAS,LTOT,
     +                                 COLMAP,
     +                                 IVEC,
     +                                 XVEC,
     +
     +                                     C (1,BASNR+1) )
     +
     +
                      END IF

                      CALL  NLO__FORM_L_SITESYM_COEFFS
     +
     +                           ( NBAS,MXLSIZE,
     +                             LTOT,LDIM,NAL,
     +                             TSYM,
     +                             P,
     +                             XVEC,
     +                             XMAT,
     +
     +                                     C (1,BASNR+1) )
     +
     +
                      CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                           ( NBAS,LTOT,
     +                             NBAS,NBAS,
     +                             NBAS,LTOT,
     +                             LTOT,
     +                             LTOT,NBAS,
     +                             0,
     +                             SAVEP,.FALSE.,.FALSE.,
     +                             C (1,BASNR+1),
     +                             P,
     +
     +                                     W (BASNR+1),
     +
     +                             XMAT )
     +
     +
                      DO J = 1,LTOT
                         WEIGHT = W (BASNR+J)
                         IF (DABS (WEIGHT).LT.WEQZERO) THEN
                             W (BASNR+J) = ZERO
                         END IF
                         BDSIZE (BASNR+J) = 1
                         BDATOM (BASNR+J) = B
                      END DO
                   END DO

                   BASL = BASL + LTOT

                END DO
            END IF
         END DO
C
C
C             ...determine the NAO angular symmetry.
C
C
         LTRG  = .TRUE.
         UTRG  = .FALSE.
         SAVEP = .TRUE.

         BASNR = 0

         DO ATOM = 1,NATOM

            NXBA = BASEND (ATOM) - BASBEG (ATOM) + 1

            CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                 ( NBAS,NXBA,
     +                   NBAS,NBAS,
     +                   NBAS,NXBA,
     +                   NXBA,
     +                   NXBA,NBAS,
     +                   0,
     +                   SAVEP,LTRG,UTRG,
     +                   C (1,BASNR+1),
     +                   P,
     +                   XVEC,
     +
     +                            XMAT )
     +
     +
            CALL  MAT__BLOCKDIAG_REAL_SYMMETRIC
     +
     +                 ( NBAS,NXBA,
     +                   MXNBA,NXBA,
     +                   NXBA,NXBA,
     +                   NXBA,
     +                   NXBA,
     +                   PSYMACC,
     +                   REVERS,
     +                   IVEC (1),
     +                   IVEC (NXBA+1),
     +                   XMAT,
     +
     +                            XVEC,
     +                            PA )
     +
     +
            NLTYPE = NSHELLS (ATOM)

            BASL = NXBA

            DO N = NLTYPE,1,-1
               NAL = NBASAL (N,ATOM)
               LTYPE = SHELLS (N,ATOM)
               LDIM = LSIZE (LTYPE)
               LTOT = NAL * LDIM

               DO J = 1,NXBA
                  X = ZERO
                  DO I = 1,LTOT
                     X = X + PA (BASL-I+1,J)**2
                  END DO
                  XVEC (J) = ZERO
                  XVEC (NXBA+J) = X
               END DO

               DO J = 1,NXBA
                  X = XVEC (NXBA+J)
                  DO I = 1,NXBA
                     XVEC (I) = XVEC (I) + X * PA (I,J)**2
                  END DO
               END DO

               DO I = 1,NXBA
                  IF (XVEC (I).GT.PSYMACC) THEN
                      ANGSYM (BASNR+I) = LTYPE
                  END IF
               END DO

               BASL = BASL - LTOT

            END DO

            BASNR = BASNR + NXBA

         END DO
C
C
C             ...printout of NAO weights and coefficient matrix for
C                testing.
C
C
C         CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
C     +
C     +                ( 6,
C     +                  ' NAO weight vector ',
C     +                  NBAS,
C     +                  NBAS,
C     +                  W )
C     +
C     +
C         CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
C     +
C     +              ( 6,
C     +                ' NAO coefficient matrix in AO basis ',
C     +                NBAS,NBAS,
C     +                NBAS,NBAS,
C     +                C )
C     +
C     +
C
C
C             ...ready!
C
C
         RETURN
         END
