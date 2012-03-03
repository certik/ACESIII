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
         SUBROUTINE  NLO__SYMMETRY_RELATED_CENTERS
     +
     +                    ( NATOM,
     +                      ATDONE,
     +                      IDXSYM,
     +                      DMAP,
     +                      SYMMAP,
     +                      DOMAT,
     +                      DSYMACC,
     +                      XYZ,
     +                      R,X,Y,Z,
     +                      DIST,
     +                      XVEC,
     +
     +                              PLATONIC,
     +                              RING,NRING,RINGSZ,
     +                              NTETRA,NCUBE,NOCTA,NDODECA,NICOSA,
     +                              PLATO )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__SYMMETRY_RELATED_CENTERS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine consists of two parts: part 1) and 2).
C
C                  Part 1): For each atomic center C in the collection
C                           of all atomic centers there exist sets
C                           of symmetry related centers which are
C                           equally distant from center C. Of all these
C                           sets, we pick out only those that have the
C                           minimum size, excluding sizes 1 and 2.
C                           Their indices will be placed into array
C                           RING in the column appropriate for center C.
C                           The # of such rings and their minimum size
C                           will be placed into the respective arrays
C                           NRING and RINGSZ at the position
C                           corresponding to C.
C
C
C                  Part 2): Determine, if the collection of all atomic
C                           centers have certain symmetry related
C                           subsets belonging to the following possible
C                           5 platonic solid symmetries:
C
C                               1) Centers arranged tetrahedrally
C                               2) Centers arranged cubically
C                               3) Centers arranged octahedrally
C                               4) Centers arranged dodecahedrally
C                               5) Centers arranged icosahedrally
C
C                           The way to do this is to analyze the
C                           'distance' info in the distance order
c                           matrix for all groups of symmetry related
C                           atomic centers.
C
C                           Procedure:
C
C                           As an example consider the centers arranged
C                           on the vertices of a cube:
C
C
C
C                                          1---------2
C                                         /|        /|
C                                        / |       / |
C                                       /  |      /  |
C                                      5---------6   |
C                                      |   |     |   |
C                                      |   3-----|---4
C                                      |  /      |  /
C                                      | /       | /
C                                      |/        |/
C                                      7---------8
C
C
C                
C                           Lets pick vertex 1 with a 'distance' value
C                           of D0. Then vertices 2,5,3 are all
C                           equidistant from vertex 1, say a distance
C                           D1. Likewise for vertices 4,6,7 by distance
C                           D2. Vertex 8 is the most distant with
C                           distance D3. Since this distance map must
C                           be obeyed by all other vertices as well,
C                           we can identify a set of 8 atomic centers
C                           belonging to the cubic symmetry if:
C
C                               i) they are all symmetry related
C                              ii) one of the centers obeys the cubic
C                                  distance map
C
C                           The above mentioned symmetry cases have
C                           then the following distance maps (D0 > D1
C                           > D2 > D3 > ...):
C
C                               1) Tetrahedron:   D0 (x1)
C                                                 D1 (x3)
C
C                               2) Cube:          D0 (x1)
C                                                 D1 (x3)
C                                                 D2 (x3)
C                                                 D3 (x1)
C
C                               3) Octahedron:    D0 (x1)
C                                                 D1 (x4)
C                                                 D2 (x1)
C
C                               4) Dodecahedron:  D0 (x1)
C                                                 D1 (x3)
C                                                 D2 (x6)
C                                                 D3 (x6)
C                                                 D4 (x3)
C                                                 D5 (x1)
C
C                               5) Icosahedron:   D0 (x1)
C                                                 D1 (x5)
C                                                 D2 (x5)
C                                                 D3 (x1)
C
C
C                Besides the distance maps, we can also specify the
C                exact values for the distance order ratios:
C
C                       [D0 - D(n+1)] / [D0 - D1]  ,  n = 0,1,2,...
C
C                The a,b distance order matrix element between two
C                atoms is defined as:
C
C                      DOMAT (a,b) = (dmax - d(ab)) / dmax
C
C                where dmax is the maximum distance in the complete
C                atomic set and d(ab) is the ordinary distance between
C                atoms a and b. From this we imediately see that all D0
C                are = 1. For the above ratio we obtain with this
C                definition, using the fact that D(n+1) denotes the
C                distance order between vertex n+1 and 0:
C
C                   [D0 - D(n+1)] / [D0 - D1]  =  d(n+1,0) / d(1,0)
C
C                a ratio involving only internal polyhedra distances.
C                These are always expressible as a linear function
C                of the polyhedra side length and hence the ratio
C                becomes independent of size and is specific for each
C                platonic solid. The following is a list of the ratios
C                expected:
C
C                    1) Tetrahedron:  (D0 - D1) / (D0 - D1) = 1
C
C                    2) Cube:         (D0 - D1) / (D0 - D1) = 1
C                                     (D0 - D2) / (D0 - D1) = sqrt (2)
C                                     (D0 - D3) / (D0 - D1) = sqrt (3)
C
C                    3) Octahedron:   (D0 - D1) / (D0 - D1) = 1
C                                     (D0 - D2) / (D0 - D1) = sqrt (2)
C
C                    4) Dodecahedron: (D0 - D1) / (D0 - D1) = 1
C                                     (D0 - D2) / (D0 - D1) = ?
C                                     (D0 - D3) / (D0 - D1) = ? 
C                                     (D0 - D4) / (D0 - D1) = ?
C                                     (D0 - D5) / (D0 - D1) = ?
C
C                    5) Icosahedron:  (D0 - D1) / (D0 - D1) = 1
C                                     (D0 - D2) / (D0 - D1) = ?
C                                     (D0 - D3) / (D0 - D1) = ?
C
C
C
C
C                  Input:
C
C                    NATOM        =  total # of atomic centers
C                    ATDONE (A)   =  will be used as an indicator if
C                                    atom A has been symmetry checked
C                                    or not.
C                    IDXSYM       =  array that will be used to store
C                                    symmetry related atomic indices.
C                    DMAP         =  array that will be used to store
C                                    distance maps, i.e. the frequency
C                                    of identical 'distances'
C                    SYMMAP (A,B) =  has been set equal to 1 if atoms
C                                    A and B were found to be symmetry
C                                    related. A value of 0 indicates
C                                    no symmetry relation.
C                    DOMAT (A,B)  =  distance order matrix containing
C                                    'distance' info between atoms A
C                                    and B.
C                    DSYMACC      =  symmetry accuracy for distance
C                                    order matrix values
C                    XYZ          =  x,y,z-coordinates for all the
C                                    atoms in 3 x NATOM matrix format.
C                    R            =  will hold atomic indices
C                    X,Y,Z        =  will hold atomic x,y,z-coordinates
C                    DIST         =  array that will be used to store
C                                    distance orders
C                    XVEC         =  flp scratch vector
C                    PLATONIC     =  is true, if the current state of
C                                    the calculation thinks the density
C                                    of the molecule possesses platonic
C                                    symmetry.
C
C
C                  Output:
C
C                    PLATONIC     =  is still true, if platonic symmetry
C                                    is assumed, false otherwise.
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

         LOGICAL     ALERT
         LOGICAL     CHECK
         LOGICAL     NSYM1,NSYM4,NSYM6,NSYM8,NSYM12,
     +               NSYM20,NSYM24,NSYM30,NSYM60
         LOGICAL     NSYMOK
         LOGICAL     PLANAR
         LOGICAL     PLATONIC
         LOGICAL     PROCEED
         LOGICAL     TETRA,CUBE,OCTA,DODECA,ICOSA
         LOGICAL     SYM4,SYM8,SYM6,SYM20,SYM12
         LOGICAL     SYMBAD1,SYMBAD2,SYMBAD3,SYMBAD4

         INTEGER     ATOM,ATOMA,ATOMB
         INTEGER     BADSIZE
         INTEGER     CENTER
         INTEGER     I,M,N
         INTEGER     MINRING,MINSIZE
         INTEGER     NATOM
         INTEGER     NBASE
         INTEGER     NDIFF
         INTEGER     NLAST
         INTEGER     NRINGS
         INTEGER     NSIZE
         INTEGER     NSYM
         INTEGER     NTETRA,NCUBE,NOCTA,NDODECA,NICOSA

         INTEGER     ATDONE  (1:NATOM)
         INTEGER     DMAP    (1:NATOM)
         INTEGER     IDXSYM  (1:NATOM)
         INTEGER     NRING   (1:NATOM)
         INTEGER     R       (1:NATOM)
         INTEGER     RINGSZ  (1:NATOM)

         INTEGER     PLATO   (1:NATOM,1:5    )
         INTEGER     RING    (1:NATOM,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM,1:NATOM)

         DOUBLE PRECISION  DIFF
         DOUBLE PRECISION  DSYMACC

         DOUBLE PRECISION  DIST (1:NATOM)
         DOUBLE PRECISION  X    (1:NATOM)
         DOUBLE PRECISION  Y    (1:NATOM)
         DOUBLE PRECISION  Z    (1:NATOM)
         DOUBLE PRECISION  XVEC (1:NATOM)

         DOUBLE PRECISION  DOMAT  (1:NATOM,1:NATOM)
         DOUBLE PRECISION  XYZ    (1:3    ,1:NATOM)
C
C
C------------------------------------------------------------------------
C
C
C             ...initialize the data.
C
C
         BADSIZE = NATOM + 1

         DO ATOMA = 1,NATOM
            ATDONE (ATOMA) = 0
            RINGSZ (ATOMA) = BADSIZE
            NRING  (ATOMA) = 0
         END DO

         CALL  MAT__C_EQ_ZERO_INTEGER
     +
     +              ( NATOM,NATOM,
     +                NATOM,NATOM,
     +
     +                         RING )
     +
     +
         TETRA  = .FALSE.
         CUBE   = .FALSE.
         OCTA   = .FALSE.
         DODECA = .FALSE.
         ICOSA  = .FALSE.

         NTETRA  = 0
         NCUBE   = 0
         NOCTA   = 0
         NDODECA = 0
         NICOSA  = 0

         CALL  MAT__C_EQ_ZERO_INTEGER
     +
     +              ( NATOM,5,
     +                NATOM,5,
     +
     +                         PLATO )
     +
     +
C
C
C             ...outer loop over all atomic centers.
C
C
         DO ATOMA = 1,NATOM

            CHECK = ATDONE (ATOMA) .EQ. 0
C
C
C             ...find total # of symmetry related atoms corresponding
C                to current checked atom.
C
C
            IF (CHECK) THEN

                NSYM = 0
                DO ATOMB = 1,NATOM
                   IF (SYMMAP (ATOMB,ATOMA) .EQ. 1) THEN
                       NSYM = NSYM + 1
                       IDXSYM (NSYM) = ATOMB
                       ATDONE (ATOMB) = 1
                   END IF
                END DO
C
C
C             ...the higher symmetry (platonic solids) cases first.
C                Pick first vertex, determine its distance order
C                vector and from it its distance order map.
C
C
                IF (PLATONIC) THEN

                    ATOMB = IDXSYM (1)

                    DO N = 1,NSYM
                       XVEC (N) = DOMAT (IDXSYM (N),ATOMB)
                    END DO

                    CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +                         ( NSYM,NSYM,
     +                           1,NSYM,
     +                           .FALSE.,.FALSE.,
     +                           0,
     +                           XVEC,
     +
     +                                   DMAP )
     +
     +
                    DO N = 1,NSYM
                       DIST (N) = XVEC (DMAP (N))
                       DMAP (N) = 0
                    END DO

                    NLAST = 1
                    NDIFF = 1

                    DO N = 2,NSYM
                       DIFF = ABS (DIST (N-1) - DIST (N))
                       IF (DIFF.GT.DSYMACC) THEN
                           NSIZE = N - NLAST
                           DMAP (NDIFF) = NSIZE
                           NDIFF = NDIFF + 1
                           NLAST = N
                       END IF
                    END DO
                    NSIZE = NSYM + 1 - NLAST
                    DMAP (NDIFF) = NSIZE
C
C
C             ...check, which platonic cases apply.
C
C
                    SYM4 =        NSYM.EQ.4
     +                      .AND. NDIFF.EQ.2
     +                      .AND. DMAP (1).EQ.1
     +                      .AND. DMAP (2).EQ.3

                    SYM8 =        NSYM.EQ.8
     +                      .AND. NDIFF.EQ.4
     +                      .AND. DMAP (1).EQ.1
     +                      .AND. DMAP (2).EQ.3
     +                      .AND. DMAP (3).EQ.3
     +                      .AND. DMAP (4).EQ.1

                    SYM6 =        NSYM.EQ.6
     +                      .AND. NDIFF.EQ.3
     +                      .AND. DMAP (1).EQ.1
     +                      .AND. DMAP (2).EQ.4
     +                      .AND. DMAP (3).EQ.1

                    SYM20 =       NSYM.EQ.20
     +                      .AND. NDIFF.EQ.6
     +                      .AND. DMAP (1).EQ.1
     +                      .AND. DMAP (2).EQ.3
     +                      .AND. DMAP (3).EQ.6
     +                      .AND. DMAP (4).EQ.6
     +                      .AND. DMAP (5).EQ.3
     +                      .AND. DMAP (6).EQ.1

                    SYM12 =       NSYM.EQ.12
     +                      .AND. NDIFF.EQ.4
     +                      .AND. DMAP (1).EQ.1
     +                      .AND. DMAP (2).EQ.5
     +                      .AND. DMAP (3).EQ.5
     +                      .AND. DMAP (4).EQ.1
C
C
C             ...add the current set of symmetry related atomic
C                centers to the corresponding column positions in
C                the PLATO array.
C
C
                    IF (SYM4) THEN

                        DO N = 1,4
                           NTETRA = NTETRA + 1
                           PLATO (NTETRA,1) = IDXSYM (N)
                        END DO

                    ELSE IF (SYM8) THEN

                        DO N = 1,8
                           NCUBE = NCUBE + 1
                           PLATO (NCUBE,2) = IDXSYM (N)
                        END DO

                    ELSE IF (SYM6) THEN

                        DO N = 1,6
                           NOCTA = NOCTA + 1
                           PLATO (NOCTA,3) = IDXSYM (N)
                        END DO

                    ELSE IF (SYM20) THEN

                        DO N = 1,20
                           NDODECA = NDODECA + 1
                           PLATO (NDODECA,4) = IDXSYM (N)
                        END DO

                    ELSE IF (SYM12) THEN

                        DO N = 1,12
                           NICOSA = NICOSA + 1
                           PLATO (NICOSA,5) = IDXSYM (N)
                        END DO

                    END IF

                    TETRA  = TETRA  .OR. SYM4
                    CUBE   = CUBE   .OR. SYM8
                    OCTA   = OCTA   .OR. SYM6
                    DODECA = DODECA .OR. SYM20
                    ICOSA  = ICOSA  .OR. SYM12

                END IF
C
C
C             ...determine the rings resulting from the current set
C                of symmetry related centers. Loop over all atomic
C                centers and determine the rings of minimum size > 2.
C
C
                DO ATOMB = 1,NATOM

                   DO N = 1,NSYM
                      XVEC (N) = DOMAT (IDXSYM (N),ATOMB)
                   END DO

                   CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +                        ( NSYM,NSYM,
     +                          1,NSYM,
     +                          .FALSE.,.FALSE.,
     +                          0,
     +                          XVEC,
     +
     +                                  DMAP )
     +
     +
                   DO N = 1,NSYM
                      DIST (N) = XVEC (DMAP (N))
                   END DO

                   NBASE = 0
                   NLAST = 1
                   MINSIZE = BADSIZE

                   DO N = 2,NSYM
                      DIFF = ABS (DIST (N-1) - DIST (N))
                      IF (DIFF.GT.DSYMACC) THEN
                          NSIZE = N - NLAST

                          IF (NSIZE.GE.3) THEN
                              IF (NSIZE.LT.MINSIZE) THEN
                                  NBASE = 0
                              END IF
                              DO M = 1,NSIZE
                                 NBASE = NBASE + 1
                                 ATOM = IDXSYM (DMAP (NLAST+M-1))
                                 R (NBASE) = ATOM
                                 X (NBASE) = XYZ (1,ATOM)
                                 Y (NBASE) = XYZ (2,ATOM)
                                 Z (NBASE) = XYZ (3,ATOM)
                              END DO
                              MINSIZE = MIN (MINSIZE,NSIZE)
                          END IF

                          NLAST = N
                      END IF
                   END DO

                   NSIZE = NSYM + 1 - NLAST

                   IF (NSIZE.GE.3) THEN
                       IF (NSIZE.LT.MINSIZE) THEN
                           NBASE = 0
                       END IF
                       DO M = 1,NSIZE
                          NBASE = NBASE + 1
                          ATOM = IDXSYM (DMAP (NLAST+M-1))
                          R (NBASE) = ATOM
                          X (NBASE) = XYZ (1,ATOM)
                          Y (NBASE) = XYZ (2,ATOM)
                          Z (NBASE) = XYZ (3,ATOM)
                       END DO
                       MINSIZE = MIN (MINSIZE,NSIZE)
                   END IF

                   NRINGS = NBASE / MINSIZE
C
C
C             ...check for proper subrings.
C
C
                   IF (NRINGS.GT.0) THEN

                       NBASE = 0

                       DO N = 1,NRINGS
                          CALL  NLO__ANALYZE_RINGS
     +
     +                               ( MINSIZE,
     +                                 R (NBASE+1),
     +                                 X (NBASE+1),
     +                                 Y (NBASE+1),
     +                                 Z (NBASE+1),
     +                                 DSYMACC,
     +                                 DMAP,
     +
     +                                       PLANAR )
     +
     +
                   WRITE (*,*) ' Ring # ',N
                   WRITE (*,*) ' Ring size ',MINSIZE
                   WRITE (*,*) ' At idx = ',(R(NBASE+I),I=1,MINSIZE)
                   WRITE (*,*) ' PLANAR = ',PLANAR

                          NBASE = NBASE + MINSIZE
                       END DO
                   END IF
C
C
C             ...minimum size of rings corresponding to current
C                atomic center within current set of symmetry
C                related centers has been found. Check, if it is
C                worth including into existing rings of current
C                atomic center.
C
C
                   MINRING = RINGSZ (ATOMB)
                   PROCEED = (NRINGS.GT.0) .AND. (MINSIZE.LE.MINRING)

                   IF (PROCEED) THEN
                       IF (MINSIZE.EQ.MINRING) THEN
                           NBASE = NRING (ATOMB)
                       ELSE
                           NBASE = 0
                           RINGSZ (ATOMB) = MINSIZE
                       END IF

                       NSIZE = NRINGS * MINSIZE

                       DO N = 1,NSIZE
                          RING (NBASE+N,ATOMB) = R (N)
                       END DO
                       NRING (ATOMB) = NBASE + NSIZE
                   END IF
C
C
C             ...next atomic center for ring check.
C
C
                END DO
C
C
C             ...next atomic center for symmetry check.
C
C
            END IF
         END DO
C
C
C             ...if everything went allright, some platonic solids
C                cannot be simultaneously present. The tetrahedron
C                can coexist with the octahedron only (this is seen
C                easily by placing a center on the middle of each
C                tetrahedron edge, of which there are six). Both
C                cube and octahedron symmetries can coexist because
C                of their dual nature. Same with the dodecahedron
C                and icosahedron pair. Below we analyze the situation
C                for incompatible coexistence of symmetries.
C                Also we analyze the dimensions of the symmetry
C                related centers. Observe that in case of platonic
C                symmetry we can only have at most one NSYM = 1
C                case and only the following NSYM dimensions:
C
C
C                         1) Tetrahedron:  NSYM = 4,4,6,12
C                         2) Cube:         NSYM = 8,6,12,24
C                         3) Octahedron:   NSYM = 6,8,12,24
C                         4) Dodecahedron: NSYM = 20,12,30,60
C                         5) Icosahedron:  NSYM = 12,20,30,60
C
C
C                The first three NSYM numbers arise (in that order) by
C                placing centers on the vertices, center of faces and
C                edge-midpoints, while the last numbers are obtained by
C                placing them anywhere on the face surfaces excluding
C                the center.
C
C                Note, that we cannot exactly pinpoint down which of
C                the 5 platonic solids we have, as the vertex numbers
C                NSYM = 4,6,8,12,20 occur multiple times. The only
C                exception is the tetrahedra, which we know for sure
C                is present if at least one NSYM is equal to 4.
C
C
C
         IF (PLATONIC) THEN

             ALERT =      (TETRA  .AND. CUBE)
     +               .OR. (TETRA  .AND. DODECA)
     +               .OR. (TETRA  .AND. ICOSA)
     +               .OR. (CUBE   .AND. DODECA)
     +               .OR. (CUBE   .AND. ICOSA)
     +               .OR. (OCTA   .AND. DODECA)
     +               .OR. (OCTA   .AND. ICOSA)

             IF (ALERT) THEN
                 WRITE (*,*) ' Symmetry inconsistency! '
                 WRITE (*,*) ' TETRA,CUBE,OCTA,DODECA,ICOSA = ',
     +                         TETRA,CUBE,OCTA,DODECA,ICOSA
                 WRITE (*,*) ' nlo__symmetry_related_centers '
                 WRITE (1,*) ' Symmetry inconsistency! '
                 WRITE (1,*) ' TETRA,CUBE,OCTA,DODECA,ICOSA = ',
     +                         TETRA,CUBE,OCTA,DODECA,ICOSA
                 WRITE (1,*) ' nlo__symmetry_related_centers '
                 STOP
             END IF

             DO ATOMA = 1,NATOM
                ATDONE (ATOMA) = 0
             END DO

             CENTER = 0

             DO ATOMA = 1,NATOM
                CHECK = ATDONE (ATOMA) .EQ. 0
                IF (CHECK) THEN
                    NSYM = 0
                    DO ATOMB = 1,NATOM
                       IF (SYMMAP (ATOMB,ATOMA) .EQ. 1) THEN
                           NSYM = NSYM + 1
                           ATDONE (ATOMB) = 1
                       END IF
                    END DO

                    IF (NSYM.EQ.1) THEN
                        CENTER = CENTER + 1
                    END IF

                    NSYM1  = NSYM .EQ. 1
                    NSYM4  = NSYM .EQ. 4
                    NSYM6  = NSYM .EQ. 6
                    NSYM8  = NSYM .EQ. 8
                    NSYM12 = NSYM .EQ. 12
                    NSYM20 = NSYM .EQ. 20
                    NSYM24 = NSYM .EQ. 24
                    NSYM30 = NSYM .EQ. 30
                    NSYM60 = NSYM .EQ. 60

                    NSYMOK =     NSYM1 .OR.NSYM4 .OR.NSYM6 .OR.NSYM8
     +                       .OR.NSYM12.OR.NSYM20.OR.NSYM24.OR.NSYM30
     +                       .OR.NSYM60

                    SYMBAD1 =  .NOT. NSYMOK
                    SYMBAD2 =  CENTER .GT. 1
                    SYMBAD3 =       (NSYM4.OR.NSYM6.OR.NSYM8)
     +                        .AND. (DODECA.OR.ICOSA)
                    SYMBAD4 =       (NSYM20.OR.NSYM30.OR.NSYM60)
     +                        .AND. (TETRA.OR.CUBE.OR.OCTA)

                    ALERT =      SYMBAD1
     +                      .OR. SYMBAD2
     +                      .OR. SYMBAD3
     +                      .OR. SYMBAD4

                    IF (ALERT) THEN
                        WRITE (*,*) ' Symmetry inconsistency! '
                        WRITE (*,*) ' SYMBAD1 = ',SYMBAD1
                        WRITE (*,*) ' SYMBAD2 = ',SYMBAD2
                        WRITE (*,*) ' SYMBAD3 = ',SYMBAD3
                        WRITE (*,*) ' SYMBAD4 = ',SYMBAD4
                        WRITE (*,*) ' nlo__symmetry_related_centers '
                        WRITE (1,*) ' Symmetry inconsistency! '
                        WRITE (1,*) ' SYMBAD1 = ',SYMBAD1
                        WRITE (1,*) ' SYMBAD2 = ',SYMBAD2
                        WRITE (1,*) ' SYMBAD3 = ',SYMBAD3
                        WRITE (1,*) ' SYMBAD4 = ',SYMBAD4
                        WRITE (1,*) ' nlo__symmetry_related_centers '
                        STOP
                    END IF

                END IF
             END DO
         END IF
C
C
C             ...everything is ok. Reset the atomic index counters
C                from individual to set counting.
C
C
         NTETRA  = NTETRA  / 4
         NCUBE   = NCUBE   / 8
         NOCTA   = NOCTA   / 6
         NDODECA = NDODECA / 20
         NICOSA  = NICOSA  / 12

         DO N = 1,NATOM
            NSIZE = RINGSZ (N)
            IF (NSIZE.EQ.BADSIZE) THEN
                RINGSZ (N) = 0
            ELSE
                NRING (N) = NRING (N) / NSIZE
            END IF
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
