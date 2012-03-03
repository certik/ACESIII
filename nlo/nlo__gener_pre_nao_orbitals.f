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
         SUBROUTINE  NLO__GENER_PRE_NAO_ORBITALS
     +
     +                    ( NBAS,NATOM,
     +                      MXSHELL,MXNAL,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      LSIZE,
     +                      BASBEG,BASEND,
     +                      SYMMAP,
     +                      ATDONE,IDXSYM,
     +                      P,S,
     +                      SAL,CAL,WAL,
     +                      SETWRYD,
     +                      WPRERYD,WRYDAT,
     +                      MJUMP,
     +
     +                              NMB,NRB,
     +                              COLMAP,
     +                              NRYDAL,
     +                              W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_PRE_NAO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms all the pre-NAO's by diagonalizing
C                the systems of equations P(AL)C(AL) = S(AL)C(AL)W(AL)
C                for all (AL) spaces. Analyze the symmetry average
C                weights W(AL) and decompose the pre-NAO space into
C                the minimal and Rydberg pre-NAO spaces. All atoms
C                related by symmetry are treated at the same time, thus
C                ensuring the same C(AL) coefficients.
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MXNAL        =  maximum size of atomic l-shell
C                                    space. The atomic l-shell space
C                                    is the total # of contractions for
C                                    an atomic l-shell.
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A (in increasing order!).
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    LSIZE (I)    =  I-th l-shell size
C                    BASBEG (A)   =  first basis index number for atom A
C                    BASEND (A)   =  last basis index number for atom A
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
C                    P            =  full NBAS x NBAS occupation matrix.
C                    S            =  full NBAS x NBAS overlap matrix.
C                    SAL,CAL,WAL  =  submatrices S(AL) and C(AL) for
C                                    overlap and density/pre-NAO coeffs
C                                    and subvector W(AL) for weights.
C                    SETWRYD      =  is true, if the Rydberg pre-NAO
C                                    weight threshold has to be set to
C                                    the initial value of WPRERYD.
C                    WPRERYD      =  the initial Rydberg pre-NAO weight
C                                    threshold. This is the initial
C                                    Rydberg weight criterion for
C                                    pre-NAO construction.
C                    WRYDAT (A)   =  the pre-NAO weight threshold below
C                                    which a pre-NAO will be considered
C                                    of Rydberg type for atom A. If
C                                    SETWRYD is true, then we reset
C                                    these values to the initial value
C                                    WPRERYD for all atoms. If SETWRYD
C                                    is false at this stage, we will
C                                    take the current values as have
C                                    been transmitted in argument.
C                    MJUMP        =  is .true., if the m values in the
C                                    m-space are ordered such that the
C                                    same m values are separated. This
C                                    keyword is necessary because some
C                                    AO basis functions are m-ordered
C                                    differently within each l-shell.
C                                    It invokes different types of
C                                    m-averaging algorithms.
C
C
C                  Output:
C
C                    NMB,NRB      =  # of minimal and Rydberg pre-NAO's
C                                    found.
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of I-th pre-NAO,
C                                    based on atomic clustering, in the
C                                    final NMB and NRB bundled pre-NAO
C                                    coefficient array C.
C                    NRYDAL (I,A) =  # of Rydberg pre-NAO's found for
C                                    the I-th atomic l-shell for atom A.
C                                    Needed later on for updating the
C                                    Rydberg pre-NAO's when generating
C                                    the final NAO's.
C                    W            =  pre-NAO weight vector.
C                    C            =  pre-NAO coefficient matrix in AO
C                                    basis.
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

         LOGICAL     ADD
         LOGICAL     LTRG
         LOGICAL     MJUMP
         LOGICAL     NEW
         LOGICAL     REVERS
         LOGICAL     SETWRYD
         LOGICAL     ZEROC

         INTEGER     A,B
         INTEGER     ATOM
         INTEGER     BASB,BASL,BASNR
         INTEGER     INDEX
         INTEGER     I,J,L,N
         INTEGER     LDIM,LTOT,LTYPE
         INTEGER     MXSHELL,MXNAL
         INTEGER     NBAS,NATOM
         INTEGER     NAL
         INTEGER     NLTYPE
         INTEGER     NMB,NRB,NRBOLD
         INTEGER     NSYM

         INTEGER     ATDONE  (1:NATOM  )
         INTEGER     BASBEG  (1:NATOM  )
         INTEGER     BASEND  (1:NATOM  )
         INTEGER     COLMAP  (1:NBAS   )
         INTEGER     IDXSYM  (1:NATOM  )
         INTEGER     LSIZE   (0:MXSHELL)
         INTEGER     NSHELLS (1:NATOM  )

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     NRYDAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM    ,1:NATOM)

         DOUBLE PRECISION  WEIGHT,WTHRSH
         DOUBLE PRECISION  WPRERYD
         DOUBLE PRECISION  X
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  W      (1:NBAS )
         DOUBLE PRECISION  WAL    (1:MXNAL)
         DOUBLE PRECISION  WRYDAT (1:NATOM)

         DOUBLE PRECISION  C   (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  S   (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  P   (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  CAL (1:MXNAL,1:MXNAL)
         DOUBLE PRECISION  SAL (1:MXNAL,1:MXNAL)

         PARAMETER    (ZERO    =  0.D0 )
         PARAMETER    (ONE     =  1.D0 )
         PARAMETER    (REVERS  = .TRUE.)
C
C
C------------------------------------------------------------------------
C
C
C             ...set initial data.
C
C
         LTRG  = .TRUE.
         ZEROC = .FALSE.

         DO A = 1,NATOM
            ATDONE (A) = 0
         END DO

         IF (SETWRYD) THEN
             DO A = 1,NATOM
                WRYDAT (A) = WPRERYD
             END DO
         END IF
C
C
C             ...outer loop over all atoms. For each atom, find its
C                group of symmetry equivalent atoms and process these
C                at the same time.
C
C
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
C
C
C
C             ...inner loop over all l-shell types for one of the
C                symmetry equivalent atoms (they must posses the same
C                angular momentum basis, otherwise they are not per
C                definition symmetry equivalent). The following steps
C                are performed for all atoms in the symmetry equivalent
C                group:
C
C                      i) for each atom calculate the m-symmetry
C                         averaged P(AL) and S(AL) matrices.
C
C                     ii) sum all P(AL) and S(AL) matrices up and
C                         form the group symmetry averaged matrices.
C
C                    iii) solve the generalized eigenvalue problem for
C                         the resulting m-symmetry and group-symmetry
C                         averaged P(AL) and S(AL) matrices.
C
C
C
                DO N = 1,NLTYPE
                   NAL = NBASAL (N,ATOM)
                   LTYPE = SHELLS (N,ATOM)
                   LDIM = LSIZE (LTYPE)
                   LTOT = NAL * LDIM

                   DO I = 1,NSYM
                      B = IDXSYM (I)
                      BASB = BASBEG (B) - 1
                      BASNR = BASB + BASL
                      ADD = I.GT.1

                      CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                           ( NBAS,LTOT,
     +                             MXNAL,MXNAL,
     +                             LTOT,LDIM,
     +                             BASNR,
     +                             MJUMP,LTRG,ADD,
     +                             P (1,BASNR+1),
     +
     +                                     CAL )
     +
     +
                      CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                           ( NBAS,LTOT,
     +                             MXNAL,MXNAL,
     +                             LTOT,LDIM,
     +                             BASNR,
     +                             MJUMP,LTRG,ADD,
     +                             S (1,BASNR+1),
     +
     +                                     SAL )
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

                   CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                        ( 6,
     +                          ' (averaged) CAL matrix (pre-NAO) ',
     +                          MXNAL,MXNAL,
     +                          NAL,0,
     +                          CAL )
     +
     +
                   CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                        ( 6,
     +                          ' (averaged) SAL matrix (pre-NAO) ',
     +                          MXNAL,MXNAL,
     +                          NAL,0,
     +                          SAL )
     +
     +
                   CALL  MAT__GEN_EIGSYS_REAL_SYMMETRIC
     +
     +                        ( MXNAL,MXNAL,
     +                          MXNAL,MXNAL,
     +                          MXNAL,
     +                          NAL,
     +                          REVERS,
     +
     +                                  WAL,
     +                                  SAL,
     +                                  CAL )
     +
     +
                   CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                        ( 6,
     +                          ' CAL eigenvalues (pre-NAO) ',
     +                          MXNAL,1,
     +                          NAL,1,
     +                          WAL )
     +
     +
                   CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                        ( 6,
     +                          ' CAL eigenvectors (pre-NAO) ',
     +                          MXNAL,MXNAL,
     +                          NAL,NAL,
     +                          CAL )
     +
     +
C
C
C             ...for the present group of symmetry equivalent atoms
C                generate the pre-NAO coefficients in the AO basis
C                together with the corresponding average weights.
C
C
                   DO I = 1,NSYM
                      B = IDXSYM (I)
                      BASB = BASBEG (B) - 1
                      BASNR = BASB + BASL

                      CALL  MAT__W_EQ_ZERO_FLOAT
     +
     +                           ( LTOT*NBAS,
     +                             LTOT*NBAS,
     +
     +                                     C (1,BASNR+1) )
     +
     +
                      CALL  NLO__FORM_M_EXPANDED_MATRIX
     +
     +                           ( MXNAL,MXNAL,
     +                             NBAS,LTOT,
     +                             NAL,LTOT,
     +                             BASNR,
     +                             MJUMP,.FALSE.,ZEROC,
     +                             CAL,
     +
     +                                     C (1,BASNR+1) )
     +
     +
                      CALL  NLO__FORM_M_EXPANDED_WEIGHTS
     +
     +                           ( MXNAL,
     +                             LTOT,
     +                             LTOT,NAL,
     +                             MJUMP,
     +                             WAL,
     +
     +                                     W (BASNR+1) )
     +
     +
                   END DO
C
C
C             ...next shell type for present symmetry equivalent
C                group of atoms.
C
C
                   BASL = BASL + LTOT

                END DO
            END IF
C
C
C             ...next symmetry equivalent group of atoms.
C
C
         END DO
C
C
C
C             ...the complete pre-NAO coefficient matrix and weight
C                vector are ready in atomic order. For all atoms
C                determine now the decomposition into minimal and
C                Rydberg spaces.
C
C                The decomposition into minimal and Rydberg spaces
C                for each atom A is governed by the critical value
C                WRYDAT (A), which has the following meaning:
C
C                          weight  >  WRYDAT (A)    minimal
C                          weight =<  WRYDAT (A)    Rydberg
C
C                As we move along, determining each atomic section of
C                the pre-NAO coefficients, we also need to keep track
C                of which column indices will belong to the NMB and
C                NRB spaces. A column map COLMAP is established, which
C                will contain the mapping of the columns from the
C                original atomic structure to the new NMB and NRB
C                clustered structure in the active sense. Pictorially,
C                this means a columnwise restructuring of the pre-NAO
C                coefficient matrix:
C
C
C                 at 1    at 2   ...             NMB          NRB
C               --------------------         ------------------------
C              |   |   |   |   | ...        |   |   |... |   |   |...
C              |   |   |   |   | ...        |   |   |... |   |   |...
C              | N | N | N | N | ...        | a | a |... | a | a |...
C              | M | R | M | R | ...    ->  | t | t |... | t | t |...
C              | B | B | B | B | ...        | 1 | 2 |... | 1 | 2 |...
C              |   |   |   |   | ...        |   |   |... |   |   |...
C              |   |   |   |   | ...        |   |   |... |   |   |...
C
C
C
C                The active column map has the following info:
C
C                  COLMAP (atomic based index) = NMB/NRB based index
C
C                Note, that while it is easy to determine the NMB
C                index, the NRB index has to be determined after the
C                we know the size of the NMB space. Hence the initial
C                COLMAP values for the NRB space must be protected
C                from identification at a later stage, which is easy
C                to do by adding the constant NBAS, since the maximum
C                value that NMB can take is NBAS.
C
C                The whole idea behind creating the column map is the
C                fact that the NMB and NRB spaces will be manipulated
C                separately, hence for efficiency reasons in handling
C                matrix operations on these spaces we will temporarily
C                bundle each kind together in later routines.
C
C
C
         NMB = 0
         NRB = 0
         BASNR = 0

         DO ATOM = 1,NATOM
            NLTYPE = NSHELLS (ATOM)
            WTHRSH = WRYDAT (ATOM)

            DO N = 1,NLTYPE
               NAL = NBASAL (N,ATOM)
               LTYPE = SHELLS (N,ATOM)
               LDIM = LSIZE (LTYPE)
               LTOT = NAL * LDIM

               NRBOLD = NRB

               DO L = 1,LTOT
                  BASNR = BASNR + 1
                  WEIGHT = W (BASNR)
                  IF (WEIGHT.GT.WTHRSH) THEN
                      NMB = NMB + 1
                      COLMAP (BASNR) = NMB
                  ELSE
                      NRB = NRB + 1
                      COLMAP (BASNR) = NRB + NBAS
                  END IF
               END DO

               NRYDAL (N,ATOM) = NRB - NRBOLD

            END DO
         END DO

C         CALL  MAT__PRINT_A_INTEGER_NOZEROS
C     +
C     +              ( 6,
C     +                ' NRYDAL matrix ',
C     +                MXSHELL+1,NATOM,
C     +                MXSHELL+1,NATOM,
C     +                NRYDAL )
C     +
C
C
C             ...remove protection value and update the NRB section of
C                the column map.
C
C
         DO I = 1,NBAS
            INDEX = COLMAP (I)
            IF (INDEX.GT.NBAS) THEN
                INDEX = INDEX - NBAS + NMB
                COLMAP (I) = INDEX
            END IF
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
