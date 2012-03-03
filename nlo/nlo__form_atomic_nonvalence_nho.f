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
         SUBROUTINE  NLO__FORM_ATOMIC_NONVALENCE_NHO
     +
     +                    ( NBAS,MXNBA,
     +                      NBOSIZE,
     +                      ATNXB,ATXIDX,ATXOFF,
     +                      NXB,NXA,
     +                      OCCNUM,BASOFF,
     +                      NAO2NHO,
     +                      P,
     +                      C,
     +                      IDX1,IDX2,
     +                      XVEC,
     +                      XMAT,
     +
     +                              BDNCEN,
     +                              BDCEN,
     +                              BDNBAS,
     +                              BDBAS,
     +                              BDOCC,
     +                              H )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_ATOMIC_NONVALENCE_NHO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine establishes the atomic nonvalence (Core
C                and Rydberg) NHOs by the following two choices:
C
C
C                Choice NAO2NHO = .true.
C                ------------------------
C                Here the NHOs are simply taken as the previously
C                found NAOs in 1 to 1 correspondence. This simply
C                means placing unit matrices into the H matrix
C                containing the expansion coeffs of the NHOs in terms
C                of the NAOs. Observe, that since the NAOs are
C                already site-symmetry adapted, the NHOs will be as
C                well.
C
C
C                Choice NAO2NHO = .false.
C                ------------------------
C                In this case, the NHOs are not simply taken as the
C                previously determined NAOs, rather the entire NAO
C                set of each atom is mixed together by diagonalizing
C                the atomic NAO occupation matrix:
C
C                    Loop over all relevant atoms having nonvalence NAOs
C                         For each such atom
C                             1) Set up the occupation matrix in
C                                NAO space
C                             2) Diagonalize that matrix
C                             3) Save the eigenvectors as the atomic
C                                nonvalence NHOs into the NHO collection
C                                matrix H and generate some additional
C                                data associated with it
C                    continue
C
C                The atomic nonvalence NHOs will also reflect the
C                correct atomic site symmetry. Before executing however
C                the diagonalization step 2), the NAO occupation matrix
C                is analyzed for matrix subblocks, which occur when
C                different symmetry components are present. Each
C                such submatrix is diagonalized separately and the
C                eigenvectors of the complete matrix are assembled
C                from those of the submatrix having zeroes elsewhere.
C                In this form the eigenvectors reflect the correct
C                symmetry, i.e. they are expanded in the proper set
C                of NAOs. Symmetry mixing does not occur when all
C                the eigenvalues of the occupation matrix are non-zero.
C                But sometimes several eigenvalues will be zero and
C                a brute force diagonalization of the entire occupation
C                matrix will mix eigenvectors belonging to different
C                symmetries.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    MXNBA        =  maximum # of NHOs per atom.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    ATNXB (A)    =  # of atomic nonvalence NAOs on
C                                    atom A having nonvalence NAOs.
C                    ATXIDX (A)   =  atomic index for atom A having
C                                    nonvalence NAOs.
C                    ATXOFF (A)   =  index offset for atomic nonvalence
C                                    NAOs for atom A having nonvalence
C                                    NAOs. This index is equal to the
C                                    total number of atomic nonvalence
C                                    NAOs on all nonvalence NAO atoms
C                                    preceeding nonvalence NAO atom A.
C                    NXB          =  total # of nonvalence NAOs.
C                    NXA          =  total # of nonvalence NAO atoms.
C                    OCCNUM       =  occupation number for present type
C                                    of atomic nonvalence NHOs. Can be
C                                    only 1 (Core) or 0 (Rydberg).
C                    BASOFF       =  global NHO basis offset index.
C                                    This number is needed to find the
C                                    global NHO basis index for the
C                                    J-th NXB-type bond (see output
C                                    below)
C                    NAO2NHO      =  is true, if the NHOs are set equal
C                                    to the previously determined NAOs.
C                                    If false, the complete atomic
C                                    NAO occupation matrix will be
C                                    diagonalized to generate the NHOs.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    C            =  section of the NAO coefficient
C                                    matrix in AO basis corresponding
C                                    to the nonvalence NAOs treated
C                                    here.
C                    IDX1,IDX2    =  integer vectors used for indexing.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C
C
C                  Output:
C
C                    BDNCEN (J)   =  # of atomic centers for the J-th
C                                    NXB-type bond. This is always = 1
C                                    here.
C                    BDCEN (1,J)  =  the atomic center index for the
C                                    J-th NXB-type bond.
C                    BDNBAS (J)   =  # of basis functions (NHOs) for
C                                    the J-th NXB-type bond. This is
C                                    always = 1 here.
C                    BDBAS (1,J)  =  the global basis (NHO) index for
C                                    the J-th NXB-type bond.
C                    BDOCC (J)    =  # of occupied levels for the J-th
C                                    NXB-type bond. Can be only 1 or 0.
C                    H (I,J)      =  MXNBA x NXB matrix containing the
C                                    atomic nonvalence NHOs. I is the
C                                    local atomic index labeling the
C                                    atomic NAOs from which the
C                                    nonvalence NHOs were constructed.
c                                    J is the global NHO index running
C                                    over all nonvalence NXB space
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

         LOGICAL     ADD
         LOGICAL     CASE1,CASE2
         LOGICAL     CHANGE
         LOGICAL     LTRG,UTRG
         LOGICAL     NAO2NHO
         LOGICAL     REVERS
         LOGICAL     SAVEP

         INTEGER     ATIDX
         INTEGER     ATOM
         INTEGER     BASOFF
         INTEGER     COL,ROW
         INTEGER     I,J,K,L,M,N
         INTEGER     LAST
         INTEGER     MXNBA
         INTEGER     NBAS
         INTEGER     NBOSIZE
         INTEGER     NSUB
         INTEGER     NSWEEP
         INTEGER     NXB,NXA
         INTEGER     NXBA
         INTEGER     OCCNUM
         INTEGER     OFF
         INTEGER     SWEEP

         INTEGER     ATNXB   (1:NXA  )
         INTEGER     ATXIDX  (1:NXA  )
         INTEGER     ATXOFF  (1:NXA  )
         INTEGER     BDNBAS  (1:NXB  )
         INTEGER     BDNCEN  (1:NXB  )
         INTEGER     BDOCC   (1:NXB  )
         INTEGER     IDX1    (1:MXNBA)
         INTEGER     IDX2    (1:MXNBA)

         INTEGER     BDBAS   (1:NBOSIZE,1:NXB)
         INTEGER     BDCEN   (1:NBOSIZE,1:NXB)

         DOUBLE PRECISION  X
         DOUBLE PRECISION  ZERO,MATZERO

         DOUBLE PRECISION  XVEC  (1:MXNBA)

         DOUBLE PRECISION  C     (1:NBAS ,1:NXB  )
         DOUBLE PRECISION  H     (1:MXNBA,1:NXB  )
         DOUBLE PRECISION  P     (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  XMAT  (1:NBAS ,1:MXNBA)

         DATA  REVERS   /.TRUE./
         DATA  MATZERO  /1.D-12/
         DATA  ZERO     /0.D0  /
C
C
C------------------------------------------------------------------------
C
C
C             ...if requested, set the NHOs equal to the NAOs.
C
C
         IF (NAO2NHO) THEN

             DO ATOM = 1,NXA
                OFF = ATXOFF (ATOM)
                NXBA = ATNXB (ATOM)
                ATIDX = ATXIDX (ATOM)

                CALL  MAT__C_EQ_UNIT_FLOAT
     +
     +                     ( MXNBA,NXBA,
     +                       NXBA,NXBA,
     +
     +                                 H (1,OFF+1) )
     +
     +
                DO I = 1,NXBA
                   BDNBAS (OFF+I) = 1
                   BDNCEN (OFF+I) = 1
                   BDCEN (1,OFF+I) = ATIDX
                   BDBAS (1,OFF+I) = BASOFF + OFF + I
                   BDOCC (OFF+I) = OCCNUM
                END DO

             END DO

             RETURN

         END IF
C
C
C             ...the complete atomic NAO occupation matrices will
C                be diagonalized. Loop over all atoms having the
C                particular nonvalence set of NAOs.
C
C
         LTRG = .TRUE.
         UTRG  = .TRUE.
         SAVEP = .TRUE.

         DO 1000 ATOM = 1,NXA

            OFF = ATXOFF (ATOM)
            NXBA = ATNXB (ATOM)
            ATIDX = ATXIDX (ATOM)

            CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                 ( NBAS,NXBA,
     +                   NBAS,NBAS,
     +                   NBAS,NXBA,
     +                   NXBA,
     +                   NXBA,NBAS,
     +                   0,
     +                   SAVEP,LTRG,UTRG,
     +                   C (1,OFF+1),
     +                   P,
     +                   XVEC,
     +
     +                           XMAT )
     +
     +
            CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                   ( 6,
     +                     ' Atomic occupation matrix ',
     +                     NBAS,NXBA,
     +                     NXBA,0,
     +                     XMAT )
     +
     +
C
C
C             ...determine the presence of submatrices. NSUB will
C                contain the # of such submatrices and IDX1 has
C                the submatrix number on those rows of the NAO
C                occupation matrix belonging together. MATZERO is
C                the threshold for when a matrix element is considered
C                to be zero.
C
C
            DO N = 1,NXBA
               IDX1 (N) = 0
            END DO

            NSUB = 0

            DO N = 1,NXBA
               IF (IDX1 (N).EQ.0) THEN

                   NSUB = NSUB + 1
                   IDX1 (N) = NSUB

                   DO I = N+1,NXBA
                      CASE1 =         IDX1 (I) .EQ. 0
                      CASE2 = ABS (XMAT (I,N)) .GT. MATZERO
                      IF (CASE1.AND.CASE2) THEN
                          IDX1 (I) = NSUB
                      END IF
                   END DO

                   CHANGE = .TRUE.
                   NSWEEP = NXBA - N

                   DO SWEEP = 1,NSWEEP
                      IF (CHANGE) THEN
                          CHANGE = .FALSE.
                          DO J = N+1,NXBA
                             ADD = .FALSE.
                             DO I = J,NXBA
                                CASE1 =         IDX1 (I) .EQ. NSUB
                                CASE2 = ABS (XMAT (I,J)) .GT. MATZERO
                                ADD = ADD .OR. (CASE1.AND.CASE2)
                             END DO
                             IF (ADD) THEN
                               DO I = J,NXBA
                                  IF (ABS (XMAT (I,J)).GT.MATZERO) THEN
                                      IDX1 (I) = NSUB
                                      CHANGE = .TRUE.
                                  END IF
                               END DO
                             END IF
                          END DO
                      END IF
                   END DO

               END IF
            END DO

            DO N = 1,NXBA
               WRITE (*,*) ' N,IDX1 = ',N,IDX1 (N)
            END DO
C
C
C             ...loop over all submatrices, assemble them, diagonalize
C                them and form the overall eigenvectors of the NAO
C                occupation matrix. If only 1 submatrix, perform
C                straightforward overall diagonalization.
C
C
            IF (NSUB.EQ.1) THEN

                CALL  MAT__C_EQ_A_FLOAT
     +
     +                     ( NBAS,NXBA,
     +                       MXNBA,NXBA,
     +                       NXBA,NXBA,
     +                       XMAT,
     +
     +                               H (1,OFF+1) )
     +
     +
                CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                     ( MXNBA,NXBA,NXBA,
     +                       NXBA,
     +                       REVERS,
     +
     +                               XVEC,
     +                               H (1,OFF+1) )
     +
     +
            ELSE

                LAST = 0

                DO K = 1,NSUB
                   M = 0
                   DO N = 1,NXBA
                      IF (IDX1 (N).EQ.K) THEN
                          M = M + 1
                          IDX2 (M) = N
                      END IF
                   END DO

                   DO J = 1,M
                      COL = IDX2 (J)
                      DO I = J,M
                         ROW = IDX2 (I)
                         H (I,OFF+LAST+J) = XMAT (ROW,COL)
                      END DO
                   END DO

                   CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                          ( 6,
     +                            ' Atomic occupation submatrix ',
     +                            MXNBA,M,
     +                            M,0,
     +                            H (1,OFF+LAST+1) )
     +
     +
                   CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                        ( MXNBA,M,M,
     +                          M,
     +                          REVERS,
     +
     +                                  XVEC (LAST+1),
     +                                  H (1,OFF+LAST+1) )
     +
     +
                   CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                          ( 6,
     +                            ' Eigenvectors ',
     +                            MXNBA,M,
     +                            M,M,
     +                            H (1,OFF+LAST+1) )
     +
     +
                   DO J = 1,M
                      L = NXBA
                      COL = OFF + LAST + J
                      DO I = M,1,-1
                         ROW = IDX2 (I)
                         DO N = L,ROW+1,-1
                            H (N,COL) = ZERO
                         END DO
                         H (ROW,COL) = H (I,COL)
                         L = ROW - 1
                      END DO
                      DO N = L,1,-1
                         H (N,COL) = ZERO
                      END DO
                   END DO

                   LAST = LAST + M

                END DO
C
C
C             ...put eigenvalues and eigenvectors of NAO occupation
C                matrix in desired order.
C
C
                DO I = 1,NXBA-1
                   K = I
                   X = XVEC (I)

                   IF (REVERS) THEN
                       DO J = I+1,NXBA
                          IF (XVEC (J).GT.X) THEN
                              K = J
                              X = XVEC (J)
                          END IF
                       END DO
                   ELSE
                       DO J = I+1,NXBA
                          IF (XVEC (J).LT.X) THEN
                              K = J
                              X = XVEC (J)
                          END IF
                       END DO
                   END IF

                   IF (K.NE.I) THEN
                       XVEC (K) = XVEC (I)
                       XVEC (I) = X
                       DO J = 1,NXBA
                          X = H (J,OFF+I)
                          H (J,OFF+I) = H (J,OFF+K)
                          H (J,OFF+K) = X
                       END DO
                   END IF

                END DO

            END IF

C            CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
C     +
C     +                   ( 6,
C     +                     ' NHO matrix Eigenvalues ',
C     +                     NXBA,
C     +                     NXBA,
C     +                     XVEC )
C     +
C     +
C            CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
C     +
C     +                   ( 6,
C     +                     ' NHO matrix Eigenvectors ',
C     +                     MXNBA,NXBA,
C     +                     NXBA,NXBA,
C     +                     H (1,OFF+1) )
C     +
C     +
            DO I = 1,NXBA
               BDNBAS (OFF+I) = 1
               BDNCEN (OFF+I) = 1
               BDCEN (1,OFF+I) = ATIDX
               BDBAS (1,OFF+I) = BASOFF + OFF + I
               BDOCC (OFF+I) = OCCNUM
            END DO

 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
