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
         SUBROUTINE  NLO__FORM_NBO_X_CENTERS
     +
     +                    ( NBAS,NHCEN,NOCC,
     +                      BONDSIZE,NBOSIZE,
     +                      MXNBA,
     +                      ATNHB,ATHOFF,
     +                      NHB,NHA,
     +                      OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,
     +                      NHOCEN,NHOIDX,
     +                      WBOND,WSTAR,
     +                      WRYD,WOCC,
     +                      C,
     +                      PH,
     +                      H,
     +                      IVEC,
     +
     +                              NBO,
     +                              NLB,NBB,NEB,NAB,NYB,
     +                              COLMAP,
     +                              B,
     +                              W,
     +                              BDSIZE,
     +                              BDATOM,
     +                              CNBO )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_NBO_X_CENTERS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine combines a specific set of NHOs to form
C                the NBOs by diagonalizing the occupation submatrix part
C                of the hybrid occupation matrix PH in terms of the
C                selected NHO set. The dimension of this submatrix will
C                be equal to the dimension NHCEN of the chosen set of
C                NHOs. Diagonalization of the submatrix will generate a
C                set of highly occupied NBOs (bond orbitals) and a set
C                of lower occupied NBOs (anti-bond orbitals). The number
C                of bond orbitals on more than one center will be
C                checked against the expected one, and an error message
C                will be issued if they don't match.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NHCEN        =  # of atomic hybrid centers to be
C                                    used for bond construction.
C                    NOCC         =  expected # of highly occupied NBOs.
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid
C                                    NAOs on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    NHB          =  total # of Hybrid NAOs
C                    NHA          =  total # of hybrid atoms
C                    OFFx         =  offsets to protect overwriting
C                                    of column map indices for the
C                                    Lone-pair, Bond, Empty-pair,
C                                    Antibond and Rydberg NBOs
C                                    (x=LP,BD,EP,AB,RY) in column
C                                    map array COLMAP.
C                    NHOCEN (I)   =  atomic center index for I-th NHO
C                                    of the set of NHCEN NHOs.
C                    NHOIDX (I)   =  global NHO index for I-th NHO
C                                    of the set of NHCEN NHOs.
C                    WBOND        =  lowest weight criterion for each
C                                    NHCEN-center bond.
C                    WSTAR        =  highest weight criterion for each
C                                    NHCEN-center antibond.
C                    WRYD         =  highest weight below which a one-
C                                    center bond is classified as an
C                                    atomic Rydberg NBO.
C                    WOCC         =  the weight limit to decide when
C                                    a NBO is considered occupied or
C                                    virtual.
C                    C            =  NBAS x NHB part of the NAO
C                                    coefficient matrix in AO basis.
C                    PH           =  NHB x NHB valence occupation matrix
C                                    over NAOs.
C                    H (I,J)      =  MXNBA x NHB matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic Hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB NHOs,
C                                    with all NHOs belonging to a
C                                    specific atomic center being
C                                    grouped together.
C                    IVEC         =  int scratch array of vector type
C
C
C                  Output:
C
C                    NBO          =  current total # of NBOs processed.
C                    NxB          =  current # of Lone-pair, Bond,
C                                    Empty-pair, Antibond and Rydberg
C                                    NBOs found (x=L,B,E,A,Y).
C                    COLMAP (I)   =  contains the current portion of the
C                                    column map in the active sense for
C                                    the NBO ordering NHB -> NLB/NBB/
C                                    NEB/NAB/NYB such that COLMAP (I)
C                                    contains the position index (with
C                                    added offsets for temporary
C                                    protection!) of the I-th NBO in
C                                    the NLB/NBB/NEB/NAB/NYB order.
C                    B (I,J)      =  NBOSIZE x NHCEN matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO.
C                    W            =  NHCEN part of the NBO weights.
C                    BDSIZE       =  full NBO bond size vector, updated
C                                    by the current NHCEN part.
C                    BDATOM       =  full NBO atomic index map array,
C                                    updated by the current NHCEN part.
C                    CNBO         =  current NHCEN part of the NBO
C                                    coefficient matrix in AO basis.
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

         LOGICAL     REVERS

         INTEGER     ATOM,ATOMI,ATOMJ
         INTEGER     BONDSIZE,NBOSIZE
         INTEGER     I,J,K,L,N
         INTEGER     KDIM,LDIM,DIM
         INTEGER     KOFF,LOFF,OFF
         INTEGER     MXNBA
         INTEGER     NAO,NBO
         INTEGER     NBAS
         INTEGER     NBBOLD
         INTEGER     NHO,NHOI,NHOJ
         INTEGER     NITER
         INTEGER     NLB,NBB,NEB,NAB,NYB
         INTEGER     NOCC
         INTEGER     NHB,NHA
         INTEGER     NHCEN
         INTEGER     OFFLP,OFFBD,OFFEP,OFFAB,OFFRY

         INTEGER     ATNHB   (1:NHA  )
         INTEGER     ATHOFF  (1:NHA  )
         INTEGER     BDSIZE  (1:NBAS )
         INTEGER     COLMAP  (1:NHB  )
         INTEGER     IVEC    (1:NHCEN)
         INTEGER     NHOCEN  (1:NHCEN)
         INTEGER     NHOIDX  (1:NHCEN)

         INTEGER     BDATOM  (1:NBAS,1:BONDSIZE)

         DOUBLE PRECISION  SUM,SUM1,SUM2
         DOUBLE PRECISION  WEIGHT,WBOND,WSTAR,WRYD,WOCC
         DOUBLE PRECISION  X,Y
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  W     (1:NHCEN )

         DOUBLE PRECISION  B     (1:NBOSIZE,1:NHCEN )
         DOUBLE PRECISION  C     (1:NBAS   ,1:NHB   )
         DOUBLE PRECISION  CNBO  (1:NBAS   ,1:NHCEN )
         DOUBLE PRECISION  H     (1:MXNBA  ,1:NHB   )
         DOUBLE PRECISION  PH    (1:NHB    ,1:NHB   )

         DATA  ONE      /1.D0/
         DATA  ZERO     /0.D0/
         DATA  REVERS   /.TRUE./
C
C
C------------------------------------------------------------------------
C
C
C             ...handle one-center case separately to obtain the
C                Lone-pairs, Empty-pairs and Rydberg NBOs.
C
C
         IF (NHCEN.EQ.1) THEN

             NHO = NHOIDX (1)
             ATOM = NHOCEN (1)
             DIM = ATNHB (ATOM)
             OFF = ATHOFF (ATOM)

             WEIGHT = ZERO
             DO L = 1,DIM
                SUM = ZERO
                DO K = 1,DIM
                   SUM = SUM + H (K,NHO) * PH (OFF+K,OFF+L)
                END DO
                WEIGHT = WEIGHT + SUM * H (L,NHO)
             END DO

             W (1) = WEIGHT
             B (1,1) = ONE

             NBO = NBO + 1
             BDSIZE (NBO) = 1
             BDATOM (NBO,1) = ATOM

             IF (WEIGHT.LT.WRYD) THEN
                 NYB = NYB + 1
                 COLMAP (NBO) = NYB + OFFRY
             ELSE IF (WEIGHT.LT.WSTAR) THEN
                 NEB = NEB + 1
                 COLMAP (NBO) = NEB + OFFEP
             ELSE IF (WEIGHT.LT.WOCC) THEN
                 WRITE (*,*) ' WARNING! High occupation Empty-pair! '
                 WRITE (*,*) ' WEIGHT = ',WEIGHT,' ( > ',WSTAR,' ) '
                 WRITE (1,*) ' WARNING! High occupation Empty-pair! '
                 WRITE (1,*) ' WEIGHT = ',WEIGHT,' ( > ',WSTAR,' ) '
                 NEB = NEB + 1
                 COLMAP (NBO) = NEB + OFFEP
             ELSE IF (WEIGHT.LT.WBOND) THEN
                 WRITE (*,*) ' WARNING! Low occupation Lone-pair! '
                 WRITE (*,*) ' WEIGHT = ',WEIGHT,' ( < ',WBOND,' ) '
                 WRITE (1,*) ' WARNING! Low occupation Lone-pair! '
                 WRITE (1,*) ' WEIGHT = ',WEIGHT,' ( < ',WBOND,' ) '
                 NLB = NLB + 1
                 COLMAP (NBO) = NLB + OFFLP
             ELSE
                 NLB = NLB + 1
                 COLMAP (NBO) = NLB + OFFLP
             END IF

             CALL  MAT__W_EQ_ZERO_FLOAT
     +
     +                  ( NBAS,
     +                    NBAS,
     +
     +                            CNBO (1,1) )
     +
     +
             DO K = 1,DIM
                NAO = OFF + K
                X = H (K,NHO)
                DO N = 1,NBAS
                   CNBO (N,1) = CNBO (N,1) + X * C (N,NAO)
                END DO
             END DO

         ELSE
C
C
C             ...the multicenter case to obtain Bonds and Antibonds.
C                Calculate lower triangle of B matrix in NHO basis
C                from the PH matrix in NAO basis. This B matrix will
C                hold temporarily the appropriate submatrix of PH.
C
C
             DO J = 1,NHCEN
                NHOJ = NHOIDX (J)
                ATOMJ = NHOCEN (J)
                LDIM = ATNHB (ATOMJ)
                LOFF = ATHOFF (ATOMJ)

                DO I = J,NHCEN
                   NHOI = NHOIDX (I)
                   ATOMI = NHOCEN (I)
                   KDIM = ATNHB (ATOMI)
                   KOFF = ATHOFF (ATOMI)

                   SUM1 = ZERO
                   DO L = 1,LDIM
                      SUM2 = ZERO
                      DO K = 1,KDIM
                         SUM2 = SUM2 + H (K,NHOI) * PH (KOFF+K,LOFF+L)
                      END DO
                      SUM1 = SUM1 + SUM2 * H (L,NHOJ)
                   END DO

                   B (I,J) = SUM1

                END DO
             END DO
C
C
C             ...diagonalize B matrix.
C
C
             CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                    ( 1,
     +                      ' B matrix ',
     +                      NBOSIZE,NHCEN,
     +                      NHCEN,0,
     +                      B )
     +
     +
             CALL  MAT__C_EQ_A_FLOAT
     +
     +                  ( NBOSIZE,NHCEN,
     +                    NBAS,NHCEN,
     +                    NHCEN,NHCEN,
     +                    B,
     +
     +                            CNBO )
     +
     +
             CALL  NLO__SYMMETRIC_JACOBI
     +
     +                  ( NBAS,NHCEN,
     +                    NBOSIZE,NHCEN,
     +                    NHCEN,
     +                    NHCEN,
     +                    1,NHCEN,
     +                    REVERS,
     +                    IVEC,
     +                    CNBO,
     +
     +                            NITER,
     +                            W,
     +                            B )
     +
     +
C             CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
C     +
C     +                  ( NBOSIZE,NHCEN,NHCEN,
C     +                    NHCEN,
C     +                    REVERS,
C     +
C     +                            W,
C     +                            B )
C     +
C     +
             CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
     +
     +                    ( 1,
     +                      ' Eigenvalues of B matrix ',
     +                      NHCEN,
     +                      NHCEN,
     +                      W )
     +
     +
             CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                    ( 1,
     +                      ' Eigenfunctions of B matrix ',
     +                      NBOSIZE,NHCEN,
     +                      NHCEN,NHCEN,
     +                      B )
     +
     +
C
C
C             ...form the current Bonds/Antibonds set.
C
C
             NBBOLD = NBB

             CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +                  ( NBAS,NHCEN,
     +                    NBAS,NHCEN,
     +
     +                            CNBO )
     +
     +
             DO J = 1,NHCEN
                NBO = NBO + 1
                BDSIZE (NBO) = NHCEN

                WEIGHT = W (J)

                IF (WEIGHT.LT.WSTAR) THEN
                    NAB = NAB + 1
                    COLMAP (NBO) = NAB + OFFAB
                ELSE IF (WEIGHT.LT.WOCC) THEN
                    WRITE (*,*) ' WARNING! High occupation Antibond! '
                    WRITE (*,*) ' WEIGHT = ',WEIGHT,' ( > ',WSTAR,' ) '
                    WRITE (1,*) ' WARNING! High occupation Antibond! '
                    WRITE (1,*) ' WEIGHT = ',WEIGHT,' ( > ',WSTAR,' ) '
                    NAB = NAB + 1
                    COLMAP (NBO) = NAB + OFFAB
                ELSE IF (WEIGHT.LT.WBOND) THEN
                    WRITE (*,*) ' WARNING! Low occupation Bond! '
                    WRITE (*,*) ' WEIGHT = ',WEIGHT,' ( < ',WBOND,' ) '
                    WRITE (1,*) ' WARNING! Low occupation Bond! '
                    WRITE (1,*) ' WEIGHT = ',WEIGHT,' ( < ',WBOND,' ) '
                    NBB = NBB + 1
                    COLMAP (NBO) = NBB + OFFBD
                ELSE
                    NBB = NBB + 1
                    COLMAP (NBO) = NBB + OFFBD
                END IF

                DO I = 1,NHCEN
                   NHO = NHOIDX (I)
                   ATOM = NHOCEN (I)
                   KDIM = ATNHB (ATOM)
                   KOFF = ATHOFF (ATOM)

                   X = B (I,J)
                   DO K = 1,KDIM
                      NAO = KOFF + K
                      Y = X * H (K,NHO)
                      DO N = 1,NBAS
                         CNBO (N,J) = CNBO (N,J) + Y * C (N,NAO)
                      END DO
                   END DO

                   BDATOM (NBO,I) = ATOM

                END DO
             END DO
C
C
C             ...check, if current # of NBO bonds matches the
C                expected results.
C
C
             IF ((NBB-NBBOLD).NE.NOCC) THEN
                 WRITE (*,*) ' # of Bond NBOs mismatch! '
                 WRITE (*,*) ' Current #, expected # = ',NBB-NBBOLD,NOCC
                 WRITE (*,*) ' nlo__form_nbo_x_centers '
                 WRITE (1,*) ' # of Bond NBOs mismatch! '
                 WRITE (1,*) ' Current #, expected # = ',NBB-NBBOLD,NOCC
                 WRITE (1,*) ' nlo__form_nbo_x_centers '
                 STOP
             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
