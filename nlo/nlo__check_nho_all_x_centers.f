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
         SUBROUTINE  NLO__CHECK_NHO_ALL_X_CENTERS
     +
     +                    ( NATOM,
     +                      BONDSIZE,NBOSIZE,
     +                      NHCEN,NHATOM,
     +                      MXNHBA,MXNBA,
     +                      ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                      RYD2HYB,
     +                      ATHIDX,INDEX,
     +                      NHB,NHA,
     +                      WBOND,WSTAR,
     +                      SAH,
     +                      PH,PHSUB,
     +                      BOMAT,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      IVEC,
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
C  OPERATION   : NLO__CHECK_NHO_ALL_X_CENTERS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine checks, if bond forming pre-NHO orbitals
C                can be constructed over NHCEN atomic hybrid centers
C                from the present hybrid occupation matrix PH. The
C                occupation matrix PH may already have been depleted of
C                several large eigenvalues corresponding to occupied
C                bond orbitals.
C
C                The procedure loops over all distinct! NHCEN atomic
C                hybrid center combinations and tries each out for
C                possible bond forming pre-NHO construction. If, at
C                any stage of this process, a failure in the bond
C                forming pre-NHO generation is issued by the calling
C                subroutines, the routine returns immediately.
C
C                Since 2-center bonds play such a fundamental role,
C                the routine considers this special case in a different
C                way. Rather than looping by brute force over all
C                atomic pairs indiscriminantly, the loop is only over
C                those atomic pairs which were found to have the
C                largest atomic bond orders. By imposing that order
C                for the 2-center bond search we eliminate all possible
C                'surprise' bonds between two centers that might pop
C                up during the brute force search, preventing 2-center
C                bonds to be formed where they actually should be. The
C                'surprise' bonds can form when the weight criterion
C                for bond formation is low and they prevent true bond
C                formation later on because they remove the bond
C                content from the occupation matrix.
C
C
C                  Input:
C
C                    NATOM        =  total # of atomic centers
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    NHCEN        =  current # of atomic hybrid centers
C                                    to be checked for NHCEN centered
C                                    bond construction.
C                    NHATOM       =  current total # of atomic hybrid
C                                    centers left for checking.
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    ATHCEN (I)   =  current hybrid atomic labels
C                                    (indices) for I-th atomic hybrid
C                                    center within the set of NHATOM
C                                    atomic hybrid centers.
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
C                    ATHIDX (I)   =  will contain NHCEN atomic labels
C                                    out of the pool of the NHATOM
C                                    atomic labels present in array
C                                    ATHCEN.
C                    INDEX (I)    =  simple index array that will hold
C                                    index numbers from 1 to NHCEN in
C                                    strictly increasing order.
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
C                    BOMAT (A,B)  =  simplified bond order matrix
C                                    containing info if atoms A and B
C                                    are considered to be bonded or not
C                    MXCHOOSE     =  maximum # of bonds selected to
C                                    be chosen. The maximum is build
C                                    from all # of chosen bonds for all
C                                    bondsizes.
C                    NCHOOSE      =  # of bonds to be chosen for bonds
C                                    of size NHCEN. Four cases:
C                                    1) = 9999 => skip search for bonds
C                                       of size NHCEN.
C                                    2) = 0 => complete search for all
C                                       possible bonds of size NHCEN
C                                       will be performed.
C                                    3) = -n => only n bonds of size
C                                       NHCEN will be searched between
C                                       those atomic indices as
C                                       provided by the CHOOSE array.
C                                    4) = +n => same as case 3) but
C                                       followed by a complete search
C                                       for all possible remaining
C                                       bonds of size NHCEN.
C                                    Priority level: 1) > 3) > 4) > 2).
C                    CHOOSE (I,N) =  contains the I-th atomic index of
C                                    the N-th chosen bond of size NHCEN.
C                                    The order of the atomic indices
C                                    is arbitrary.
C                    IVEC         =  int scratch array of vector type
C                    XVEC         =  flp scratch array of vector type
C                    XMAT         =  flp scratch array of matrix type
C                              
C
C
C                  Output:
C
C                    FAILED       =  is true, if any problems were
C                                    found such that the present
C                                    and subordinate routines cannot
C                                    proceed.
C                    NBOND        =  current total # of bonds.
C                    NHYB (A)     =  current total # of pre-NHO's on
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
         LOGICAL     ACCEPT
         LOGICAL     CHOOSED
         LOGICAL     FAILED
         LOGICAL     INCRESE
         LOGICAL     NLO__ACCEPT_ATOMS_BOND_SEARCH
         LOGICAL     RYD2HYB

         INTEGER     BONDSIZE,NBOSIZE
         INTEGER     I,K,N
         INTEGER     IXPOS
         INTEGER     MXCHOOSE
         INTEGER     MXNHBA,MXNBA
         INTEGER     NATOM
         INTEGER     NBOND
         INTEGER     NCHOOSE
         INTEGER     NHCEN
         INTEGER     NHATOM
         INTEGER     NHB,NHA

         INTEGER     ATNHB   (1:NHA   )
         INTEGER     ATHCEN  (1:NHATOM)
         INTEGER     ATHIDX  (1:NHCEN )
         INTEGER     ATHOFF  (1:NHA   )
         INTEGER     ATHVAL  (1:NHA   )
         INTEGER     BDNBAS  (1:NHB   )
         INTEGER     BDNCEN  (1:NHB   )
         INTEGER     BDOCC   (1:NHB   )
         INTEGER     INDEX   (1:NHCEN )
         INTEGER     IVEC    (1:2*NHB+NHA)
         INTEGER     NHYB    (1:NHA   )

         INTEGER     BDBAS   (1:NBOSIZE,1:NHB)
         INTEGER     BDCEN   (1:NBOSIZE,1:NHB)
         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE)

         DOUBLE PRECISION  WBOND,WSTAR

         DOUBLE PRECISION  W    (1:NHB)
         DOUBLE PRECISION  XVEC (1:2*NHB)

         DOUBLE PRECISION  BOMAT (1:NATOM ,1:NATOM )
         DOUBLE PRECISION  H     (1:MXNBA ,1:NHB   )
         DOUBLE PRECISION  PH    (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHDEP (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHSUB (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  SAH   (1:MXNHBA,1:MXNHBA)
         DOUBLE PRECISION  XMAT  (1:NHB   ,1:NHB   )
C
C
C------------------------------------------------------------------------
C
C
C             ...if bonds of size NHCEN should be skipped, return.
C
C
         IF (NCHOOSE.EQ.9999) THEN
             RETURN
         END IF
C
C
C             ...if bonds of size NHCEN have been choosen, do these
C                first.
C
C
         CHOOSED = .FALSE.

         N = IABS (NCHOOSE)

         IF (N.GT.0) THEN

             ABSOLUT = .FALSE.
             INCRESE = .TRUE.

             DO K = 1,N

                ACCEPT = NLO__ACCEPT_ATOMS_BOND_SEARCH
     +
     +                        ( NATOM,
     +                          BONDSIZE,
     +                          NHCEN,NHA,
     +                          CHOOSE (1,K),ATHVAL,
     +                          NHYB,
     +                          BOMAT,
     +                          CHOOSED,
     +                          MXCHOOSE,NCHOOSE,CHOOSE,
     +                          IVEC )
     +
     +
                IF (ACCEPT) THEN

                    CALL  NLO__SORT_INT_VECTOR_ELEMENTS
     +
     +                         ( NHCEN,NHCEN,
     +                           1,NHCEN,
     +                           ABSOLUT,INCRESE,
     +                           0,
     +                           CHOOSE (1,K),
     +
     +                                    IVEC )
     +
     +
                    DO I = 1,NHCEN
                       ATHIDX (I) = CHOOSE (IVEC (I),K)
                    END DO

                    CALL  NLO__CHECK_NHO_X_CENTERS
     +
     +                         ( NHCEN,
     +                           NBOSIZE,
     +                           MXNHBA,MXNBA,
     +                           ATHIDX,ATNHB,ATHVAL,ATHOFF,
     +                           RYD2HYB,
     +                           NHB,NHA,
     +                           WBOND,WSTAR,
     +                           SAH,
     +                           PH,PHSUB,
     +                           IVEC (1),IVEC (NHB+1),IVEC (2*NHB+1),
     +                           XVEC (1),
     +                           XVEC (NHB+1),XMAT,
     +
     +                                    FAILED,
     +                                    NBOND,
     +                                    NHYB,
     +                                    BDNCEN,
     +                                    BDCEN,
     +                                    BDNBAS,
     +                                    BDBAS,
     +                                    BDOCC,
     +                                    W,H,
     +                                    PHDEP )
     +
     +
                    IF (FAILED) THEN
                        RETURN
                    END IF
                END IF

             END DO

             IF (NCHOOSE.LT.0) THEN
                 RETURN
             END IF

             CHOOSED = .TRUE.

         END IF
C
C
C             ...the complete search for all bonds of size NHCEN.
C                Handle special situation for 1-center bonds.
C
C
         IF (NHCEN.EQ.1) THEN

             DO K = 1,NHATOM

                ATHIDX (1) = ATHCEN (K)

                CALL  NLO__CHECK_NHO_X_CENTERS
     +
     +                     ( 1,
     +                       NBOSIZE,
     +                       MXNHBA,MXNBA,
     +                       ATHIDX,ATNHB,ATHVAL,ATHOFF,
     +                       RYD2HYB,
     +                       NHB,NHA,
     +                       WBOND,WSTAR,
     +                       SAH,
     +                       PH,PHSUB,
     +                       IVEC (1),IVEC (NHB+1),IVEC (2*NHB+1),
     +                       XVEC (1),
     +                       XVEC (NHB+1),XMAT,
     +
     +                                FAILED,
     +                                NBOND,
     +                                NHYB,
     +                                BDNCEN,
     +                                BDCEN,
     +                                BDNBAS,
     +                                BDBAS,
     +                                BDOCC,
     +                                W,H,
     +                                PHDEP )
     +
     +
                IF (FAILED) THEN
                    RETURN
                END IF

             END DO

             RETURN

         END IF
C
C
C             ...the general case for bond with >= 2 centers. Set
C                initial atomic hybrid center index combination:
C                1 2 3 4 5 ...
C
C
         ABSOLUT = .FALSE.
         INCRESE = .TRUE.

         DO I = 1,NHCEN
            IVEC (I) = ATHCEN (I)
            INDEX (I) = I
         END DO

         ACCEPT = NLO__ACCEPT_ATOMS_BOND_SEARCH
     +
     +                 ( NATOM,
     +                   BONDSIZE,
     +                   NHCEN,NHA,
     +                   IVEC,ATHVAL,
     +                   NHYB,
     +                   BOMAT,
     +                   CHOOSED,
     +                   MXCHOOSE,NCHOOSE,CHOOSE,
     +                   IVEC (NHCEN+1) )
     +
     +
         IF (ACCEPT) THEN

             CALL  NLO__SORT_INT_VECTOR_ELEMENTS
     +
     +                  ( NHCEN,NHCEN,
     +                    1,NHCEN,
     +                    ABSOLUT,INCRESE,
     +                    0,
     +                    IVEC,
     +
     +                             IVEC (NHCEN+1) )
     +
     +
             DO I = 1,NHCEN
                ATHIDX (I) = IVEC (IVEC (NHCEN+I))
             END DO

             CALL  NLO__CHECK_NHO_X_CENTERS
     +
     +                  ( NHCEN,
     +                    NBOSIZE,
     +                    MXNHBA,MXNBA,
     +                    ATHIDX,ATNHB,ATHVAL,ATHOFF,
     +                    RYD2HYB,
     +                    NHB,NHA,
     +                    WBOND,WSTAR,
     +                    SAH,
     +                    PH,PHSUB,
     +                    IVEC (1),IVEC (NHB+1),IVEC (2*NHB+1),
     +                    XVEC (1),
     +                    XVEC (NHB+1),XMAT,
     +
     +                             FAILED,
     +                             NBOND,
     +                             NHYB,
     +                             BDNCEN,
     +                             BDCEN,
     +                             BDNBAS,
     +                             BDBAS,
     +                             BDOCC,
     +                             W,H,
     +                             PHDEP )
     +
     +
             IF (FAILED) THEN
                 RETURN
             END IF
         END IF

         K = 0
C
C
C             ...find all combinations of the NHCEN index numbers
C                such that each index number is greater than the
C                preceeding one. If we encounter an index combination
C                that gives us a hybrid atom combination involving
C                at least one atom with an already complete bond forming
C                pre-NHO set, that index combination is not accepted
C                for x-center bond search, and we move to the next
C                possible one. The same, if the index combination
C                indicates isolated (not bonded) atoms, as judged by
C                the simplified bond order matrix.
C
C
 1000    IXPOS = NHCEN - K
         INDEX (IXPOS) = INDEX (IXPOS) + 1

         IF (INDEX (IXPOS) .LE. (NHATOM-K)) THEN

             DO I = 1,K
                INDEX (IXPOS+I) = INDEX (IXPOS) + I
             END DO

             DO I = 1,NHCEN
                IVEC (I) = ATHCEN (INDEX (I))
             END DO

             ACCEPT = NLO__ACCEPT_ATOMS_BOND_SEARCH
     +
     +                     ( NATOM,
     +                       BONDSIZE,
     +                       NHCEN,NHA,
     +                       IVEC,ATHVAL,
     +                       NHYB,
     +                       BOMAT,
     +                       CHOOSED,
     +                       MXCHOOSE,NCHOOSE,CHOOSE,
     +                       IVEC (NHCEN+1) )
     +
     +
C
C
C             ...if the current atomic index pattern is found to
C                be accepted, reorder the atomic indices such that
C                indices i,j,k,... will be ordered in ascending
C                order i<j<k<...to allow easy construction of the
C                corresponding density submatrix. Note, that the
C                atomic index vector ATHIDX is being overwritten
C                by this reordering.
C
C
             IF (ACCEPT) THEN

                 CALL  NLO__SORT_INT_VECTOR_ELEMENTS
     +
     +                      ( NHCEN,NHCEN,
     +                        1,NHCEN,
     +                        ABSOLUT,INCRESE,
     +                        0,
     +                        IVEC,
     +
     +                                 IVEC (NHCEN+1) )
     +
     +
                 DO I = 1,NHCEN
                    ATHIDX (I) = IVEC (IVEC (NHCEN+I))
                 END DO

                 CALL  NLO__CHECK_NHO_X_CENTERS
     +
     +                      ( NHCEN,
     +                        NBOSIZE,
     +                        MXNHBA,MXNBA,
     +                        ATHIDX,ATNHB,ATHVAL,ATHOFF,
     +                        RYD2HYB,
     +                        NHB,NHA,
     +                        WBOND,WSTAR,
     +                        SAH,
     +                        PH,PHSUB,
     +                        IVEC (1),IVEC (NHB+1),IVEC (2*NHB+1),
     +                        XVEC (1),
     +                        XVEC (NHB+1),XMAT,
     +
     +                                 FAILED,
     +                                 NBOND,
     +                                 NHYB,
     +                                 BDNCEN,
     +                                 BDCEN,
     +                                 BDNBAS,
     +                                 BDBAS,
     +                                 BDOCC,
     +                                 W,H,
     +                                 PHDEP )
     +
     +
                 IF (FAILED) THEN
                     RETURN
                 END IF
             END IF

             K = 0
         ELSE
             K = K + 1
             IF (K .GT. (NHCEN-1)) THEN
                 RETURN
             END IF
         END IF

         GOTO 1000
C
C
C             ...ready!
C
C
         END 
