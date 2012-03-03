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
         SUBROUTINE  NLO__INITIALIZE_RUN
     +
     +                    ( NBAS,NATOM,
     +                      MXSHELL,
     +                      MAXOCC,BONDSIZE,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      SPHERIC,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      INDEX,
     +                      XYZ,ZATOM,
     +                      S,
     +                      XVEC,
     +                      XMAT,
     +
     +                            NO2CEN,
     +                            WPRERYD,WNAOVAL,WNAORYD,WNBOOCC,
     +                            WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                            DSYMACC,LSYMACC,PSYMACC,QSYMACC,
     +                            MXNBA,MXLSIZE,
     +                            LSIZE,
     +                            BASBEG,BASEND,
     +                            SYMMAP,
     +                            BDATOM,
     +                            GPOINT,
     +                            DOMAT,
     +                            SOMAT,
     +                            BOMAT,
     +                            PLATONIC,
     +                            D )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__INITIALIZE_RUN
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine initializes the present run. It generates
C                all data, which is needed several times later on.
C                Also transmitted data is checked for consistency and
C                the program is aborted, if any meaningless data is
C                found. Initialize appropriate arrays to zero.
C
C                This routine generates the distance order matrix, the
C                atomic overlap order matrix, the atomic bond order
C                matrix and the occupation matrix from the density
C                matrix (i.e. the expansion coefficients of the density
C                rho(r) in terms of basis function products) and the
C                overlap matrix. The occupation matrix is defined as a
C                matrix representation of the density rho(r) in terms
C                of the basis functions.
C
C                Procedure:
C
C                The definition of the density matrix with elements
C                D (i,j) in terms of the basis functions is:
C
C                       rho (r)  =  sum  D (i,j) * chi (i) * chi (j)
C                                    ij
C
C                The occupation matrix element P (k,l) is thus:
C
C                       P (k,l)  =  < chi (k) | rho (r) | chi (l) >
C
C                                =  sum  S (k,i) * D (i,j) * S (j,l)
C                                    ij
C
C                where S is the overlap matrix of the basis functions.
C
C                           ------------------------------
C
C                The atomic bond order matrix is defined as summations
C                of elements of D*S between atomic subblocks a and b:
C
C                    BOMAT (a,b) =  sum     sum    DS (i,j) * DS (j,i)
C                                 i in a  j in b
C
C                Note that the product D*S is needed for evaluation
C                of both matrices, hence the first step is the
C                generation of D*S. The atomic bond order matrix can
C                be used to decide initially if platonic symmetry is
C                present in the density matrix. Simply diagonalize
C                the atomic bond order matrix and check for the
C                resulting eigenvalues for degeneracies > 2.
C
C                           ------------------------------
C
C                The atomic overlap order matrix is finally defined as
C                summations of square elements of S between atomic
C                subblocks a and b:
C
C                    SOMAT (a,b) =  sum     sum    S (i,j) * S (i,j)
C                                 i in a  j in b
C
C                This matrix is further normalized such that all
C                diagonal elements become equal to 1 via:
C
C                                        _______________
C                     SOMAT' (a,b) =   \/ (1/SOMAT(a,a))
C                                        _______________
C                                    * \/ (1/SOMAT(b,b))
C
C                                    * SOMAT (a,b)
C
C
C                           ------------------------------
C
C                The distance order matrix is a matrix reflecting
C                the distances between atoms in the following way:
C
C                      i) The atom itself has a distance order
C                         of 1.
C
C                     ii) The distance order between those atoms
C                         which are most appart is set to 0.
C
C                This is achieved by defining the distance order
C                between atoms A and B as:
C
C                      DOMAT (A,B) = (dmax - d(AB)) / dmax
C
C                where dmax is the maximum distance in the complete
C                atomic set and d(AB) is the ordinary distance between
C                atoms A and B. The importance of the distance order
C                matrix lies in the fact that similar arrangements
C                of atoms differing only by all distances scaled by
C                a specific factor (i.e. imploded or exploded shapes
C                from the original geometry) lead to the same distance
C                order matrix. Also note that the distance orders are
C                independent of units used (Angstrom or Bohr).
C
C                           ------------------------------
C
C                Since the density matrix D is not needed later on,
C                the array D will be overwritten with the occupation
C                matrix.
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A (in increasing order!).
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    SPHERIC      =  is true, if the l-shells are
C                                    spherical, false if they are
C                                    cartesian.
C                    MXCHOOSE     =  maximum # of bonds selected to
C                                    be chosen. The maximum is build
C                                    from all # of chosen bonds for all
C                                    bondsizes.
C                    NCHOOSE (B)  =  # of bonds to be chosen for bonds
C                                    of size B. Four cases:
C                                    1) = 9999 => skip search for bonds
C                                       of size B.
C                                    2) = 0 => complete search for all
C                                       possible bonds of size B will
C                                       be performed.
C                                    3) = -n => only n bonds of size B
C                                       will be searched between those
C                                       atomic indices as provided by
C                                       the CHOOSE array.
C                                    4) = +n => same as case 3) but
C                                       followed by a complete search
C                                       for all possible remaining
C                                       bonds of size B.
C                                    Priority level: 1) > 3) > 4) > 2).
C                    CHOOSE       =  element CHOOSE (I,N,B) contains
C                                    the I-th atomic index of the N-th
C                                    chosen bond of size B. The order
C                                    of the atomic indices is arbitrary.
C                    INDEX        =  will hold info about used atoms
C                    XYZ          =  x,y,z-coordinates for all the
C                                    atoms in 3 x NATOM matrix format.
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    S            =  normalized overlap matrix
C                    XVEC         =  flp scratch vector
C                    XMAT         =  flp scratch matrix
C                    D            =  normalized density matrix
C
C                  Output:
C
C                    NO2CEN       =  2-center bond formation criterion
C                                    for analysis of the atomic bond
C                                    order matrix.
C                    WPRERYD      =  the initial Rydberg pre-NAO weight
C                                    threshold. This is the initial
C                                    Rydberg weight criterion for
C                                    pre-NAO construction.
C                    WNAOVAL      =  the Valence NAO weight threshold.
C                    WNAORYD      =  the Rydberg NAO weight threshold.
C                    WNBOOCC      =  the weight limit to decide when
C                                    a NBO is considered occupied or
C                                    virtual.
C                    WBDMIN       =  initial minimum accepted weight
C                                    for bond construction.
C                    WSTMAX       =  initial maximum accepted weight
C                                    for antibond construction.
C                    WBDCRT       =  critical weight limit for bond
C                                    construction.
C                    WSTEP        =  weight stepping size to decrease/
C                                    increase the bond/antibond weight
C                                    limits.
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
C                    MXNBA        =  maximum # of AOs per atom
C                    MXLSIZE      =  maximum l-shell size, i.e. largest
C                                    m-subspace.
C                    LSIZE (I)    =  I-th l-shell size
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
C                    SYMMAP (A,B) =  this matrix will contain the info,
C                                    if atoms A and B are symmetry
C                                    related (=1) or not (=0). Initially
C                                    we know that for sure A is symmetry
C                                    related to A, hence we initiallize
C                                    this matrix with a unit matrix.
C                    BDATOM       =  Natural localized orbital atomic
C                                    index map array. Will be set = 0.
C                    GPOINT       =  if the center of mass of the
C                                    atomic arrangement coincides with
C                                    one of the atoms, this variable
C                                    will contain the corresponding
C                                    atomic index. If no atom coincides
C                                    with the center of mass location
C                                    a value of 0 is transmitted.
C                    DOMAT        =  distance order matrix
C                    SOMAT        =  atomic overlap order matrix
C                    BOMAT        =  atomic bond order matrix
C                    PLATONIC     =  is true, if the diagonalization
C                                    of the atomic bond order matrix
C                                    gave eigenvalues of degeneracies
C                                    > 2, showing that the density
C                                    matrix might correspond to a
C                                    platonic symmetry (ruling out
C                                    accidential degeneracies!)
C                    D            =  occupation matrix in AO basis
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

         LOGICAL     PLATONIC
         LOGICAL     SPHERIC

         INTEGER     A,B,I,J,N
         INTEGER     ATOM
         INTEGER     BASNR
         INTEGER     BONDSIZE
         INTEGER     GPOINT
         INTEGER     IBEG,IEND,JBEG,JEND
         INTEGER     LTYPE
         INTEGER     MAXDEG
         INTEGER     MXCHOOSE
         INTEGER     MXLSIZE
         INTEGER     MXNBA
         INTEGER     MXSHELL
         INTEGER     NBAS,NATOM
         INTEGER     NDEG
         INTEGER     NLAST
         INTEGER     NSHELL
         INTEGER     NXBA

         INTEGER     BASBEG  (1:NATOM)
         INTEGER     BASEND  (1:NATOM)
         INTEGER     INDEX   (1:NATOM)
         INTEGER     LSIZE   (0:MXSHELL)
         INTEGER     NCHOOSE (1:BONDSIZE)
         INTEGER     NSHELLS (1:NATOM)
         INTEGER     ZATOM   (1:NATOM)

         INTEGER     BDATOM  (1:NBAS     ,1:BONDSIZE)
         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM   )
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM   )
         INTEGER     SYMMAP  (1:NATOM    ,1:NATOM   )

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE,1:BONDSIZE)

         DOUBLE PRECISION  DAB,DMAX
         DOUBLE PRECISION  DIFF
         DOUBLE PRECISION  DSYMACC,LSYMACC,PSYMACC,QSYMACC
         DOUBLE PRECISION  MASS,MOLMASS
         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  NO2CEN
         DOUBLE PRECISION  WBDMIN,WSTMAX,WBDCRT,WSTEP
         DOUBLE PRECISION  WPRERYD,WNAORYD,WNAOVAL,WNBOOCC
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  ZERO,ONE,SMALL,THRESH

         DOUBLE PRECISION  XVEC  (1:NATOM)

         DOUBLE PRECISION  BOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  D     (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  DOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  S     (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  SOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  XMAT  (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  XYZ   (1:3    ,1:NATOM)

         PARAMETER    (ZERO    = 0.D0 )
         PARAMETER    (ONE     = 1.D0 )
         PARAMETER    (SMALL   = 1.D-6)
         PARAMETER    (THRESH  = 1.D-8)
C
C
C------------------------------------------------------------------------
C
C
C             ...set the weight limits and the symmetry accuracy
C                parameters.
C
C
         NO2CEN  = 0.025D0 * MAXOCC
         WPRERYD = 0.25D0  * MAXOCC
         WNAORYD = 0.050D0 * MAXOCC
         WNAOVAL = 0.99D0  * MAXOCC
         WNBOOCC = 0.500D0 * MAXOCC
         WBDMIN  = 0.975D0 * MAXOCC
         WSTMAX  = 0.05D0  * MAXOCC
         WBDCRT  = 0.75D0  * MAXOCC
         WSTEP   = 0.01D0  * MAXOCC

         DSYMACC = 1.D-10
         LSYMACC = 1.D-8
         PSYMACC = 1.D-10
         QSYMACC = 1.D-10
C
C
C             ...predetermine l-shell sizes needed and set maximum.
C
C
         LSIZE (0) = 1
         IF (SPHERIC) THEN
             DO LTYPE = 1,MXSHELL
                LSIZE (LTYPE) = LSIZE (LTYPE-1) + 2
             END DO
         ELSE
             WRITE (*,*) ' Cannot run with cartesian AO basis! '
             WRITE (*,*) ' SPHERIC = ',SPHERIC
             WRITE (*,*) ' nlo__initialize_run '
             WRITE (1,*) ' Cannot run with cartesian AO basis! '
             WRITE (1,*) ' SPHERIC = ',SPHERIC
             WRITE (1,*) ' nlo__initialize_run '
             STOP
         END IF

         MXLSIZE = LSIZE (MXSHELL)
C
C
C             ...determine atomic basis indices and maximum # of
C                AOs on for further use. Check here, if the l-shell
C                types are in increasing order, i.e. s,p,d,f,...
C
C
         MXNBA = 0
         BASNR = 0

         DO A = 1,NATOM
            BASBEG (A) = BASNR + 1
            NSHELL = NSHELLS (A)
            LTYPE = SHELLS (1,A)
            BASNR = BASNR + NBASAL (1,A) * LSIZE (LTYPE)

            DO N = 2,NSHELL
               IF (SHELLS (N,A).GE.LTYPE) THEN
                   LTYPE = SHELLS (N,A)
                   BASNR = BASNR + NBASAL (N,A) * LSIZE (LTYPE)
               ELSE
                   WRITE (*,*) ' L-shells not in increasing order! '
                   WRITE (*,*) ' For atom # ',A,' there is: '
                   WRITE (*,*) ' Shells # ',N-1,' is = ',LTYPE
                   WRITE (*,*) ' Shells # ',N,' is = ',SHELLS (N,A)
                   WRITE (*,*) ' nlo__initialize_run '
                   WRITE (1,*) ' L-shells not in increasing order! '
                   WRITE (1,*) ' For atom # ',A,' there is: '
                   WRITE (1,*) ' Shells # ',N-1,' is = ',LTYPE
                   WRITE (1,*) ' Shells # ',N,' is = ',SHELLS (N,A)
                   WRITE (1,*) ' nlo__initialize_run '
                   STOP
               END IF
            END DO

            BASEND (A) = BASNR
            NXBA = BASEND (A) - BASBEG (A) + 1
            MXNBA = MAX (MXNBA,NXBA)
         END DO
C
C
C             ...check the selected atomic indices for bond formation
C                (if any). Stop, if any inconsistencies are found.
C
C
         DO B = 1,BONDSIZE

            N = IABS (NCHOOSE (B))

            IF (N.NE.9999) THEN

                IF (N.GT.MXCHOOSE) THEN
                    WRITE (*,*) ' Dim for CHOOSE array too small! '
                    WRITE (*,*) ' B,N,MXCHOOSE = ',B,N,MXCHOOSE
                    WRITE (*,*) ' nlo__initialize_run '
                    WRITE (1,*) ' Dim for CHOOSE array too small! '
                    WRITE (1,*) ' B,N,MXCHOOSE = ',B,N,MXCHOOSE
                    WRITE (1,*) ' nlo__initialize_run '
                    STOP
                ELSE IF (N.GT.0) THEN

                    DO I = 1,N
                       DO ATOM = 1,NATOM
                          INDEX (ATOM) = 0
                       END DO
                       DO A = 1,B
                          ATOM = CHOOSE (A,I,B)
                          INDEX (ATOM) = INDEX (ATOM) + 1
                       END DO
                       DO ATOM = 1,NATOM
                          IF (INDEX (ATOM).GT.1) THEN
                              WRITE (*,*) ' Atom idx used > once! '
                              WRITE (*,*) ' ATOM,I,B = ',ATOM,I,B
                              WRITE (*,*) ' nlo__initialize_run '
                              WRITE (1,*) ' Atom idx used > once! '
                              WRITE (1,*) ' ATOM,I,B = ',ATOM,I,B
                              WRITE (1,*) ' nlo__initialize_run '
                              STOP
                          END IF
                       END DO
                    END DO

                END IF
            END IF
         END DO
C
C
C             ...set the initial atomic symmetry map matrix equal to
C                a unit matrix.
C
C
         DO B = 1,NATOM
            DO A = 1,NATOM
               SYMMAP (A,B) = 0
            END DO
            SYMMAP (B,B) = 1
         END DO
C
C
C             ...set natural localized orbital atomic index array
C                to zero.
C
C
         CALL  MAT__C_EQ_ZERO_INTEGER
     +
     +              ( NBAS,BONDSIZE,
     +                NBAS,BONDSIZE,
     +
     +                        BDATOM )
     +
     +
C
C
C             ...find out, if one of the atom locations coincides
C                with the center of mass location.
C
C
         XB = ZERO
         YB = ZERO
         ZB = ZERO
         MOLMASS = ZERO

         DO A = 1,NATOM
            XA = XYZ (1,A)
            YA = XYZ (2,A)
            ZA = XYZ (3,A)
            MASS = DFLOAT (ZATOM (A))
            XB = XB + MASS * XA
            YB = YB + MASS * YA
            ZB = ZB + MASS * ZA
            MOLMASS = MOLMASS + MASS
         END DO

         XB = XB / MOLMASS
         YB = YB / MOLMASS
         ZB = ZB / MOLMASS

         N = 0
         GPOINT = 0

         DO A = 1,NATOM
            XA = XYZ (1,A)
            YA = XYZ (2,A)
            ZA = XYZ (3,A)
            XC = XA - XB
            YC = YA - YB
            ZC = ZA - ZB
            DAB = DSQRT (XC*XC + YC*YC + ZC*ZC)
            IF (DAB.LT.SMALL) THEN
                N = N + 1
                IF (N.GT.1) THEN
                    WRITE (*,*) ' 2 atoms at center of mass! '
                    WRITE (*,*) ' Atom indices = ',GPOINT,A
                    WRITE (*,*) ' nlo__initialize_run '
                    WRITE (1,*) ' 2 atoms at center of mass! '
                    WRITE (1,*) ' Atom indices = ',GPOINT,A
                    WRITE (1,*) ' nlo__initialize_run '
                    STOP
                ELSE
                    GPOINT = A
                END IF
            END IF
         END DO
C
C
C             ...form D*S next.
C
C
         CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,NBAS,
     +                D,S,
     +
     +                         XMAT )
     +
     +
C
C
C             ...evaluate the atomic bond order matrix.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NATOM,NATOM,
     +                NATOM,NATOM,
     +
     +                         BOMAT )
     +
     +
         DO B = 1,NATOM
            JBEG = BASBEG (B)
            JEND = BASEND (B)
            DO A = 1,NATOM
               IBEG = BASBEG (A)
               IEND = BASEND (A)

               DO I = IBEG,IEND
               DO J = JBEG,JEND
                  BOMAT (A,B) = BOMAT (A,B) + XMAT (I,J) * XMAT (J,I)
               END DO
               END DO

            END DO
         END DO
C
C
C             ...form occupation matrix P = S * DS overwriting D.
C
C
         CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,NBAS,
     +                S,XMAT,
     +
     +                         D )
     +
     +
C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  'printing SPS with P the HF',
C     +                  NBAS,NBAS,
C     +                  NBAS,NBAS,
C     +                  D )
C     +
C     +
C
C
C             ...diagonalize the bond order matrix and determine
C                the maximum degeneracy of its eigenvalues. Use
C                the still undetermined DOMAT array for scratch.
C
C
         CALL    MAT__C_EQ_A_FLOAT
     +
     +                ( NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  BOMAT,
     +
     +                          DOMAT )
     +
     +
         CALL    MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                ( NATOM,NATOM,NATOM,
     +                  NATOM,
     +                  .TRUE.,
     +
     +                          XVEC,
     +                          DOMAT )
     +
     +
         NLAST = 1
         MAXDEG = 0

         DO N = 2,NATOM
            DIFF = ABS (XVEC (N-1) - XVEC (N))
            IF (DIFF.GT.DSYMACC) THEN
                NDEG = N - NLAST
                MAXDEG = MAX (MAXDEG,NDEG)
                NLAST = N
            END IF
         END DO
         NDEG = NATOM + 1 - NLAST
         MAXDEG = MAX (MAXDEG,NDEG)

         PLATONIC = MAXDEG .GT. 2
C
C
C             ...evaluate the distance order matrix.
C
C
         DMAX = ZERO

         DO B = 1,NATOM
            XB = XYZ (1,B)
            YB = XYZ (2,B)
            ZB = XYZ (3,B)
            DO A = B+1,NATOM
               XA = XYZ (1,A)
               YA = XYZ (2,A)
               ZA = XYZ (3,A)
               XC = XA - XB
               YC = YA - YB
               ZC = ZA - ZB
               DAB = DSQRT (XC*XC + YC*YC + ZC*ZC)
               DMAX = MAX (DAB,DMAX)
               DOMAT (A,B) = DAB
            END DO
         END DO

         DO B = 1,NATOM
            DOMAT (B,B) = ONE
            DO A = B+1,NATOM
               DOMAT (A,B) = (DMAX - DOMAT (A,B)) / DMAX
               DOMAT (B,A) = DOMAT (A,B)
            END DO
         END DO
C
C
C             ...evaluate the atomic overlap order matrix.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NATOM,NATOM,
     +                NATOM,NATOM,
     +
     +                         SOMAT )
     +
     +
         DO B = 1,NATOM
            JBEG = BASBEG (B)
            JEND = BASEND (B)
            DO A = 1,NATOM
               IBEG = BASBEG (A)
               IEND = BASEND (A)

               DO I = IBEG,IEND
               DO J = JBEG,JEND
                  SOMAT (A,B) = SOMAT (A,B) + S (I,J) * S (I,J)
               END DO
               END DO

            END DO
         END DO

         DO A = 1,NATOM
            SOMAT (A,A) = ONE / DSQRT (SOMAT (A,A))
         END DO

         DO B = 1,NATOM
            XB = SOMAT (B,B)
            DO A = B+1,NATOM
               SOMAT (A,B) = SOMAT (A,A) * XB * SOMAT (A,B)
               SOMAT (B,A) = SOMAT (A,B)
            END DO
            SOMAT (B,B) = ONE
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
