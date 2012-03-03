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
         SUBROUTINE  NLO__CHECK_NHO_X_CENTERS
     +
     +                    ( NHCEN,
     +                      NBOSIZE,
     +                      MXNHBA,MXNBA,
     +                      ATHIDX,ATNHB,ATHVAL,ATHOFF,
     +                      RYD2HYB,
     +                      NHB,NHA,
     +                      WBOND,WSTAR,
     +                      SAH,
     +                      PH,PHSUB,
     +                      INDEX,USED,DEPLET,
     +                      OVLAVG,
     +                      XVEC,XMAT,
     +
     +                              FAILED,
     +                              NBOND,
     +                              NHYB,
     +                              BDNCEN,
     +                              BDCEN,
     +                              BDNBAS,
     +                              BDBAS,
     +                              BDOCC,
     +                              W,H,
     +                              PHDEP )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__CHECK_NHO_X_CENTERS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine checks, if bond forming pre-NHO orbitals
C                can be constructed from the present hybrid occupation
C                matrix PH spreading over fixed NHCEN atomic hybrid
C                centers.
C
C                Procedure (rough sketch, more details in code below):
C
C                    i) Set up submatrix PHSUB of PH by taking
C                       appropriate subblocks of PH matrix
C                       corresponding to the fixed NHCEN atomic
C                       hybrid centers combination.
C
C                   ii) Diagonalize PHSUB and consider those
C                       eigenfunctions which have eigenvalues
C                       > WBOND.
C
C                  iii) For each such eigenfunction find the NHCEN-1
C                       partner eigenfunction as those having maximum
C                       average atomic overlap and check if they
C                       correspond to eigenvalues > WBOND or < WSTAR.
C                       If such a partner set can be found for
C                       the original eigenfunction, accept this
C                       eigenfunction for atomic pre-NHO decomposition
C                       (splitting of eigenfunction into the individual
C                       normalized Hybrid NAO components for each atom
C                       to define the pre-NHOs).
C
C                   iv) Deplete PHSUB from all eigenvalues > WBOND
C                       corresponding to the eigenfunctions + their
C                       partners found in steps ii) and iii) used
C                       for NHO decomposition.
C
C                    v) Replace PHSUB section of PH by the depleted
C                       PHSUB matrix, if depletion was performed at all.
C
C
C                  Input:
C
C                    NHCEN        =  # of atomic hybrid centers to be
C                                    checked for bond construction.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    ATHIDX (I)   =  atomic label (index) for I-th
C                                    atomic hybrid center out of the
C                                    NHCEN atomic centers entering a
C                                    NHCEN centered bond construction.
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid NAOs
C                                    on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    should be included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                                    If false, only the Valence NAO
C                                    space will be used.
C                    NHB          =  total # of Hybrid NAOs.
C                    NHA          =  total # of hybrid atoms.
C                    WBOND        =  lowest weight above which a
C                                    NHCEN-centered bond formation
C                                    is accepted.
C                    WSTAR        =  highest weight below which a
C                                    NHCEN-centered antibond formation
C                                    is accepted.
C                    SAH          =  will contain the atomic hybrid
C                                    overlap matrices between pre-NHOs.
C                    PH           =  current (possibly frequently large
C                                    eigenvalue depleted) NHB x NHB
C                                    hybrid occupation matrix.
C                    PHSUB        =  will contain the submatrix of
C                                    the occupation matrix PH
C                                    corresponding to a specific
C                                    NHCEN atomic hybrid center
C                                    combination.
C                    INDEX (I)    =  index for ordering the average
C                                    atomic overlap values.
C                    USED (I)     =  indicator for which eigenfunctions
C                                    have been used during the PHSUB
C                                    eigenfunctions analysis.
C                    DEPLET (I)   =  index for which eigenfunctions
C                                    of the PHSUB matrix have to be
C                                    used for depletion from the PHSUB
C                                    matrix.
C                    OVLAVG (I)   =  will contain the average atomic
C                                    overlap values.
C                    XVEC         =  flp scratch array of vector type
C                    XMAT         =  flp scratch array of matrix type
C
C
C                  Output:
C
C                    FAILED       =  is true, if any problems were
C                                    found such that the present
C                                    and subordinate routines cannot
C                                    proceed.
C                    NBOND        =  current total # of bonds.
C                    NHYB (A)     =  current total # of pre-NHOs on
C                                    each hybrid atom A.
C                    BDNCEN (J)   =  # of atomic centers for J-th bond.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond.
C                    BDNBAS (J)   =  # of basis functions (NHOs) for
C                                    J-th bond.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond.
C                    BDOCC (J)    =  # of occupied levels for J-th bond.
C                    W            =  pre-NHO weight vector to accumulate
C                                    the weights for WSW procedure.
C                    H (I,J)      =  MXNBA x NHB matrix containing the
C                                    current found set of pre-NHOs for
C                                    each hybrid atom. I is the local
C                                    atomic index labeling the atomic
C                                    hybrid NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB space
C                                    pre-NHOs in that order, with all
C                                    pre-NHOs belonging to a specific
C                                    atomic center being grouped
C                                    together.
C                    PHDEP        =  used for accumulation of the 
C                                    depleted occupation matrix PH.
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
         LOGICAL     LTRG
         LOGICAL     REVERS
         LOGICAL     RYD2HYB

         INTEGER     ATOM
         INTEGER     I,J,K,L,M,N
         INTEGER     MLARGE
         INTEGER     MXNHBA,MXNBA
         INTEGER     NBOND
         INTEGER     NBOSIZE
         INTEGER     NCOUNT
         INTEGER     NDEP,NPAR
         INTEGER     NHOIDX
         INTEGER     NHYBA
         INTEGER     NHB,NHBA,NHA
         INTEGER     NHCEN
         INTEGER     NITER
         INTEGER     NVAL
         INTEGER     OFF
         INTEGER     ROWX

         INTEGER     ATHIDX  (1:NHCEN)
         INTEGER     ATNHB   (1:NHA  )
         INTEGER     ATHOFF  (1:NHA  )
         INTEGER     ATHVAL  (1:NHA  )
         INTEGER     BDNBAS  (1:NHB  )
         INTEGER     BDNCEN  (1:NHB  )
         INTEGER     BDOCC   (1:NHB  )
         INTEGER     DEPLET  (1:NHCEN)
         INTEGER     INDEX   (1:NHB  )
         INTEGER     NHYB    (1:NHA  )
         INTEGER     USED    (1:NHB  )

         INTEGER     BDBAS   (1:NBOSIZE,1:NHB)
         INTEGER     BDCEN   (1:NBOSIZE,1:NHB)

         DOUBLE PRECISION  DENOM
         DOUBLE PRECISION  E,F
         DOUBLE PRECISION  EFRAC
         DOUBLE PRECISION  NORM
         DOUBLE PRECISION  OVLSUM
         DOUBLE PRECISION  SJJ,SKJ,SKK
         DOUBLE PRECISION  VALSQR,VALTOL
         DOUBLE PRECISION  WBOND,WSTAR,WHIGH
         DOUBLE PRECISION  XJ,XK
         DOUBLE PRECISION  ZERO,ONE,TWO,TINY

         DOUBLE PRECISION  OVLAVG (1:NHB)
         DOUBLE PRECISION  W      (1:NHB)
         DOUBLE PRECISION  XVEC   (1:NHB)

         DOUBLE PRECISION  H      (1:MXNBA ,1:NHB   )
         DOUBLE PRECISION  PH     (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHDEP  (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHSUB  (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  SAH    (1:MXNHBA,1:MXNHBA)
         DOUBLE PRECISION  XMAT   (1:NHB   ,1:NHB   )

         DATA  ONE      /1.D0/
         DATA  TWO      /2.D0/
         DATA  ZERO     /0.D0/
         DATA  EFRAC    /1.D-6/
         DATA  TINY     /1.D-10/
         DATA  VALTOL   /0.5D0/
         DATA  REVERS   /.TRUE./
C
C
C------------------------------------------------------------------------
C
C
C             ...extract lower triangle of PHSUB matrix from PH matrix.
C
C
         WRITE (1,*) ' Atomic Centers = ',(ATHIDX(I),I=1,NHCEN)

         LTRG = .TRUE.
         FAILED = .FALSE.
         EXTRACT = .TRUE.

         CALL  NLO__HANDLE_MATRIX_SECTIONS
     +
     +              ( NHB,NHB,
     +                NHB,NHB,
     +                NHA,NHCEN,
     +                NHCEN,
     +                ATHIDX,ATNHB,ATHOFF,
     +                EXTRACT,
     +                LTRG,
     +                PH,
     +
     +                        M,
     +                        PHSUB )
     +
     +
C         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
C     +
C     +                ( 1,
C     +                  ' PHSUB matrix ',
C     +                  NHB,NHB,
C     +                  M,0,
C     +                  PHSUB )
C     +
C     +
C
C
C             ...if the Rydberg space has been added to the Valence
C                space for the pre-NHO construction and the number
C                of centers > 2, there is a danger of accidental
C                degeneracies present in the PHSUB matrix, which
C                will result in a mixing of eigenvectors belonging
C                to different symmetry species. The presence of the
C                Rydberg space can lower the eigenvalues of the
C                highest occupied states to almost exactly = 2 even
C                if they belong to different symmetries. A way out
C                of this is to remove a tiny fraction of an electron
C                on the diagonal elements which are highest in
C                absolute value for each atom. That way one introduces
C                a negligible negative shift on these eigenvalues
C                removing the accidental degeneracy.
C
C                The idea is based on including a controlled
C                perturbation on the eigenvalues of PHSUB. The key
C                result is Eq.(44.10) of Wilkinsons book "The
C                Algebraic Eigenvalue Problem", Oxford 1965, which
C                states that when a matrix B (not necessarily small)
C                is added to matrix A as a perturbation, then the
C                eigenvalues of A are changed by an amount which
C                lies between the smallest and greatest of the
C                eigenvalues of B. In our case we know the exact
C                eigenvalues of the pertrubation matrix B, since
C                we define it to be diagonal with zeroes or very
C                small negative numbers -X. Hence if the original
C                eigenvalues of PHSUB are denoted by p(i) and the
C                perturbed ones by p'(i), we then have:
C
C
C                        p(i)-|X(max)| =< p'(i) =< p(i)
C
C
C                The fraction of electron to be removed on each
C                atomic center should be large enough to remove
C                the degeneracy but small enough such that it will
C                not affect further searches for atomic NHOs. Since
C                each atomic center can have only a very limited
C                number of bond forming NHOs, applying the perturbation
C                a couple of times to the overall PH matrix will
C                result in a removal of several tiny fractions of
C                electrons from the PH matrix at each atomic center.
C                An electron fraction of 1.D-6 is an appropriate
C                value for this purpose.
C
C
C
         IF (RYD2HYB .AND. NHCEN.GT.2) THEN

             DO N = 1,NHCEN
                ATOM = ATHIDX (N)
                INDEX (N) = ATNHB (ATOM)
             END DO

             CALL  NLO__PERTURB_OCCUPATION_MATRIX
     +
     +                  ( NHB,NHB,
     +                    M,NHCEN,
     +                    INDEX,
     +                    EFRAC,
     +                    USED,XVEC,
     +
     +                             PHSUB )
     +
     +
         END IF
C
C
C             ...diagonalize PHSUB matrix.
C
C
         CALL  NLO__SYMMETRIC_JACOBI
     +
     +              ( NHB,NHB,
     +                NHB,NHB,
     +                NHB,
     +                M,
     +                1,M,
     +                REVERS,
     +                USED,
     +                PHSUB,
     +
     +                         NITER,
     +                         XVEC,
     +                         XMAT )
     +
     +
C         CALL  MAT__C_EQ_A_FLOAT
C     +
C     +              ( NHB,NHB,
C     +                NHB,NHB,
C     +                M,M,
C     +                PHSUB,
C     +
C     +                        XMAT )
C     +
C     +
C         CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
C     +
C     +              ( NHB,NHB,NHB,
C     +                M,
C     +                REVERS,
C     +
C     +                        XVEC,
C     +                        XMAT )
C     +
C     +
C
C
C             ...we now identify those eigenfunctions, which have
C                eigenvalues > WBOND. Let us take the first such
C                eigenfunction Psi (1). Since it is formaly an
C                expansion over one NHO at each atomic center,
C                we can write:
C
C                        Psi (1) = sum  c (i) NHO (i)
C                                   i
C
C                where i runs over all NHCEN atomic centers. Psi (1)
C                has now (NHCEN - 1) partner functions:
C
C
C                        Psi' (x) = sum  c' (i) NHO' (i)
C                                    i
C
C                which however do not have the exact structure of
C                the NHOs as for Psi (1), that is:
C
C                        NHO' (i) = only approximately NHO (i)
C
C                The goal is thus to identify which eigenfunctions
C                of the whole set obey this NHO approximation to
C                a maximum extent. For this, we look at the NHCEN
C                atomic overlaps of each normalized (to eliminate the
C                c(i) and c'(i) dependence!) atomic NHO and NHO'
C                for any pair of eigenfunctions Psi (1) and Psi (x),
C                where x runs over all the eigenfunctions different
C                from Psi (1) and determine the average over all the
C                atomic overlaps. As a result we will obtain a set of
C                x average atomic overlaps for each Psi (1)/Psi (x)
C                pair. The (NHCEN - 1) eigenfunctions Psi (x)
C                corresponding to the largest average atomic overlaps
C                are then identified as the (NHCEN - 1) partner
C                functions of Psi (1).
C
C                All partners of Psi (1) correspond to (NHCEN - 1)
C                eigenvalues, which, for Psi (1) to be accepted
C                for NHO decomposition, should either be > WBOND
C                or < WSTAR. If any eigenvalue of the partners is found
C                to ly outside that range, we discard Psi (1) for NHO
C                decomposition and look at the next remaining unused
C                Psi (2) in the set of eigenfunctions with eigenvalues
C                > WBOND.
C                
C                First, select the high weight space of eigenfunctions
C                with eigenvalues > WBOND.
C
C
C
         MLARGE = 0
         DO J = 1,M
            IF (XVEC (J) .GT. WBOND) THEN
                MLARGE = MLARGE + 1
            END IF
         END DO
C
C
C             ...before we start checking the high weight eigenfunctions
C                for their partners in case of NHCEN > 1, we must be
C                sure that we order them in such a way that those
C                which will be checked for NHO decomposition have
C                actually all c(i) expansion coefficients nonzero.
C                Otherwise we might run into the problem of
C                generating a zero pre-NHO on an atomic site, which
C                is unacceptable. Since we start checking the
C                eigenfunctions for NHO decomposition from the first
C                one on, it is sufficient to order (if necessary)
C                them in such a way that those high weight
C                eigenfunctions with all their c(i) nonzero come
C                first. This is checked by forming the sum of squares
C                of the atomic valence part expansion coefficients
C                for each eigenfunction and moving the latter if
C                necessary.
C
C
         IF (NHCEN.GT.1) THEN

             DO J = 1,MLARGE
                ROWX = 0
                OVLSUM = ONE
                DO N = 1,NHCEN
                   ATOM = ATHIDX (N)
                   NHBA = ATNHB (ATOM)
                   NVAL = ATHVAL (ATOM)

                   SJJ = ZERO
                   DO I = 1,NVAL
                      XJ = XMAT (ROWX+I,J)
                      SJJ = SJJ + XJ * XJ
                   END DO
                   OVLSUM = DMIN1 (OVLSUM,SJJ)
                   ROWX = ROWX + NHBA
                END DO
                OVLAVG (J) = OVLSUM
             END DO

             ABSOLUT = .FALSE.
             INCRESE = .FALSE.

             CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +                  ( MLARGE,MLARGE,
     +                    1,MLARGE,
     +                    ABSOLUT,INCRESE,
     +                    0,
     +                    OVLAVG,
     +
     +                             USED )
     +
     +
             DO J = 1,MLARGE
                INDEX (USED (J)) = J
             END DO

             CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +                  ( NHB,NHB,
     +                    NHB,
     +                    NHB,
     +                    NHB,
     +                    M,MLARGE,
     +                    INDEX,
     +                    USED,
     +                    OVLAVG,
     +
     +                             XMAT )
     +
     +
             CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +                  ( NHB,
     +                    NHB,
     +                    NHB,
     +                    MLARGE,
     +                    INDEX,
     +                    OVLAVG,
     +
     +                             XVEC )
     +
     +
             CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
     +
     +                    ( 1,
     +                      ' Eigenvalues of PHSUB matrix ',
     +                      NHB,
     +                      M,
     +                      XVEC )
     +
     +
C             CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                    ( 1,
C     +                      ' Eigenfunctions of PHSUB matrix ',
C     +                      NHB,NHB,
C     +                      M,M,
C     +                      XMAT )
C     +
C     +
         END IF
C
C
C             ...loop over the high weight space eigenfunctions.
C
C
         DO J = 1,MLARGE
            USED (J) = 0
         END DO

         ABSOLUT = .FALSE.
         INCRESE = .FALSE.

         DO 1000 J = 1,MLARGE

            WHIGH = XVEC (J)
            IF (USED (J).EQ.0) THEN

                NDEP = 1
                DEPLET (1) = J
C
C
C             ...if the number of atomic centers is > 1, determine
C                all average valence atomic overlaps between the
C                current eigenfunction J and all the others.
C                Note, that if an atomic overlap is found to be
C                very small, this is an indication that one of the
C                corresponding c's in the expansion:
C
C
C                            Psi = sum  c (i) NHO (i)
C                                   i
C
C                is nearly equal to zero. In such a case we do
C                not include the corresponding atomic overlap in
C                the determination of the average.
C
C                The atomic overlaps will be determined for the
C                atomic Valence space only, as the Rydberg space
C                (if present) has only complimentary character
C                (i.e. small expansion coefficients) for the
C                partner functions of interest. Hence an eigenfunction
C                will be rejected for atomic overlap determination,
C                if its Valence expansion coefficient squares sum
C                up to less than a Valence tolerance value VALTOL.
C                These kind of eigenfunctions are dominated by the
C                Rydberg space and not of interest for partner search.
C
C
                IF (NHCEN.GT.1) THEN

                    DO K = 1,M
                       OVLAVG (K) = ZERO
                       IF (K.NE.J) THEN

                           ROWX = 0
                           VALSQR = ZERO
                           DO N = 1,NHCEN
                              ATOM = ATHIDX (N)
                              NHBA = ATNHB (ATOM)
                              NVAL = ATHVAL (ATOM)
                              DO I = 1,NVAL
                                 XK = XMAT (ROWX+I,K)
                                 VALSQR = VALSQR + XK * XK
                              END DO
                              ROWX = ROWX + NHBA
                           END DO

                           IF (VALSQR.GT.VALTOL) THEN
                               ROWX = 0
                               NCOUNT = 0
                               OVLSUM = ZERO
                               DO N = 1,NHCEN
                                  ATOM = ATHIDX (N)
                                  NHBA = ATNHB (ATOM)
                                  NVAL = ATHVAL (ATOM)

                                  SJJ = ZERO
                                  SKJ = ZERO
                                  SKK = ZERO
                                  DO I = 1,NVAL
                                     XJ = XMAT (ROWX+I,J)
                                     XK = XMAT (ROWX+I,K)
                                     SJJ = SJJ + XJ * XJ
                                     SKJ = SKJ + XK * XJ
                                     SKK = SKK + XK * XK
                                  END DO

                                  DENOM = SJJ * SKK
                                  IF (DENOM.GT.TINY) THEN
                                      NCOUNT = NCOUNT + 1
                                      OVLSUM = OVLSUM + (SKJ*SKJ)/DENOM
                                  END IF

                                  ROWX = ROWX + NHBA
                               END DO

                               IF (NCOUNT.GT.0) THEN
                                   OVLSUM = OVLSUM / DFLOAT (NCOUNT)
                               END IF

                               OVLAVG (K) = DMIN1 (OVLSUM+TINY,ONE)
                           END IF
                       END IF
                    END DO
C
C
C             ...sort all M average atomic overlaps in decreasing
C                order. Use the INDEX array to store the ordered
C                eigenfunction indices. Select the (NHCEN - 1) highest
C                average atomic overlaps as partners to the J-th
C                eigenfunction and test their eigenvalues. If an
C                eigenvalue passes the test, mark the corresponding
C                eigenfunction as used. If also the eigenvalue
C                is > WBOND mark the eigenfunction as ready for
C                depletion.
C
C
                    CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +                         ( M,M,
     +                           1,M,
     +                           ABSOLUT,INCRESE,
     +                           0,
     +                           OVLAVG,
     +
     +                                   INDEX )
     +
     +
                    NPAR = NHCEN - 1

                    DO N = 1,NPAR
                       K = INDEX (N)
                       E = XVEC (K)
                       IF (E.GT.WBOND .OR. DABS (E).LT.WSTAR) THEN
                           USED (K) = 1
                           IF (K.LE.MLARGE) THEN
                               NDEP = NDEP + 1
                               DEPLET (NDEP) = K
                           END IF
                       ELSE
                           GOTO 1000
                       END IF
                    END DO

                END IF
C
C
C             ...we get here, if the correct number of partners
C                with the desired properties have been found for
C                the current J-th eigenfunction. Deplete the J-th
C                eigenvalue and all its partners with eigenvalues
C                > WBOND from the PHSUB matrix.
C
C
                EXTRACT = .TRUE.

                CALL  NLO__HANDLE_MATRIX_SECTIONS
     +
     +                     ( NHB,NHB,
     +                       NHB,NHB,
     +                       NHA,NHCEN,
     +                       NHCEN,
     +                       ATHIDX,ATNHB,ATHOFF,
     +                       EXTRACT,
     +                       LTRG,
     +                       PHDEP,
     +
     +                               M,
     +                               PHSUB )
     +
     +
                DO N = 1,NDEP
                   I = DEPLET (N)
                   E = XVEC (I)
                   DO L = 1,M
                      F = E * XMAT (L,I)
                      DO K = L,M
                         PHSUB (K,L) = PHSUB (K,L) - F * XMAT (K,I)
                      END DO
                   END DO
                END DO

                EXTRACT = .FALSE.

                CALL  NLO__HANDLE_MATRIX_SECTIONS
     +
     +                     ( NHB,NHB,
     +                       NHB,NHB,
     +                       NHA,NHCEN,
     +                       NHCEN,
     +                       ATHIDX,ATNHB,ATHOFF,
     +                       EXTRACT,
     +                       LTRG,
     +
     +                               PHDEP,
     +
     +                       M,
     +                       PHSUB )
     +
     +
C
C
C             ...incorporate the present atomic parts of the J-th
C                averaged eigenfunction as an additional bond forming
C                pre-NHO on each atom, including their common weight.
C                If the # of bond forming pre-NHOs on any hybrid atom
C                exceeds the # of valence NAOs on that atom, the
C                routine is considered to have failed and we return.
C
C

                NBOND = NBOND + 1
                BDOCC (NBOND) = NDEP
                BDNBAS (NBOND) = NHCEN
                BDNCEN (NBOND) = NHCEN

                ROWX = 0
                DO 2000 N = 1,NHCEN
                   ATOM = ATHIDX (N)
                   OFF = ATHOFF (ATOM)
                   NVAL = ATHVAL (ATOM)
                   NHBA = ATNHB (ATOM)
                   NHYBA = NHYB (ATOM)

                   NHYBA = NHYBA + 1

                   IF (NHYBA .GT. NVAL) THEN
                       WRITE (*,*) ' # of bond pre-NHO > # of val NAO! '
                       WRITE (*,*) ' At#,#NHO,#NAO = ',ATOM,NHYBA,NVAL
                       WRITE (*,*) ' nlo__check_nho_x_centers '
                       WRITE (1,*) ' # of bond pre-NHO > # of val NAO! '
                       WRITE (1,*) ' At#,#NHO,#NAO = ',ATOM,NHYBA,NVAL
                       WRITE (1,*) ' nlo__check_nho_x_centers '
                       FAILED = .TRUE.
                       RETURN
                   END IF

                   NHOIDX = OFF + NHYBA

                   NORM = ZERO
                   DO I = 1,NHBA
                      NORM = NORM + XMAT (ROWX+I,J)**2
                   END DO
                   NORM = ONE / DSQRT (NORM)

                   DO I = 1,NHBA
                      H (I,NHOIDX) = NORM * XMAT (ROWX+I,J)
                   END DO

                   W (NHOIDX) = WHIGH
C
C
C             ...we now have to check, if the current (possibly still
C                incomplete) set of bond forming pre-NHOs on the
C                present hybrid atom is still a linearly independent
c                set. If not, we have failed and return. If everything
C                is ok, save the pre-NHO data for further use.
C
C
                   CALL  NLO__CHECK_LINEAR_DEPENDENCY
     +
     +                        ( MXNBA,NHYBA,
     +                          MXNHBA,MXNHBA,
     +                          NHBA,NHYBA,
     +                          H (1,OFF+1),
     +                          SAH,
     +
     +                                 DEPEND )
     +
     +
                   IF (DEPEND) THEN
                       WRITE (*,*) ' Linear dep. bond pre-NHO set! '
                       WRITE (*,*) ' Atom number = ',ATOM
                       WRITE (*,*) ' nlo__check_nho_x_centers '
                       WRITE (1,*) ' Linear dep. bond pre-NHO set! '
                       WRITE (1,*) ' Atom number = ',ATOM
                       WRITE (1,*) ' nlo__check_nho_x_centers '
                       FAILED = .TRUE.
                       RETURN
                   END IF

                   BDCEN (N,NBOND) = ATOM
                   BDBAS (N,NBOND) = NHOIDX

                   NHYB (ATOM) = NHYBA
                   ROWX = ROWX + NHBA

 2000           CONTINUE

            END IF
C
C
C             ...check next high weight space eigenfunction.
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
