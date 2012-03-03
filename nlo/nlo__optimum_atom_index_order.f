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
         SUBROUTINE  NLO__OPTIMUM_ATOM_INDEX_ORDER
     +
     +                    ( NATOM,N2CEN,
     +                      INDEX,
     +                      AT2CEN,
     +
     +                             ATORD )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__OPTIMUM_ATOM_INDEX_ORDER
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : Given a set of ordered atomic pair (2 center) indices,
C                this routine extracts an optimum atomic index sequence
C                via the following algorithm:
C
C                    1) Take first index pair I,J and start the
C                       atomic index sequence with I.
C
C                    2) Go through the index pair list from the
C                       beginning to the end, adding all atomic
C                       indices {K} in that order, that contain
C                       atom I in all index pairs {I,K}.
C
C                    3) Take second element I of already established
C                       partial atomic index sequence and repeat
C                       step 2). Then take third element I, etc...
C
C                Observe, that isolated sequences of atoms (i.e. not
C                bonded) might occur in cases the program deals
C                with a molecule composed of well separated fragments.
C                Also if the set of all index pairs does not contain
C                all atom indices, the atomic index sequence will be
C                completed by those missing.
C
C                A special case arises, if no index pairs were found
C                from the bond order matrix. This indicates a isolated
C                atom or a set of well separated atoms. In both these
C                cases we have N2CEN = 0 and no info is sitting in
C                array AT2CEN. Thus the case N2CEN = 0 has to be dealt
C                with separately.
C
C                  Input:
C
C                    NATOM        =  total # of atoms
C                    N2CEN        =  # of ordered atomic index pairs
C                                    that will be considered for
C                                    optimum atomic index sequence
C                                    construction.
C                    INDEX        =  will hold info about used atoms
C                    AT2CEN (1,N) =  1st atomic index of N-th ordered
C                                    pair.
C                    AT2CEN (2,N) =  2nd atomic index of N-th ordered
C                                    pair.
C
C
C                  Output:
C
C                    ATORD (I)    =  I-th atomic index in sequence.
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

         LOGICAL     CASE1,CASE2
         LOGICAL     CHECK
         LOGICAL     IGNORE

         INTEGER     AT1,AT2
         INTEGER     ATNEW
         INTEGER     ATOM
         INTEGER     I,J,M
         INTEGER     IDX
         INTEGER     N2CEN
         INTEGER     NATOM
         INTEGER     USED,NOTUSED

         INTEGER     ATORD (1:NATOM)
         INTEGER     INDEX (1:NATOM)

         INTEGER     AT2CEN (1:2,1:N2CEN)
C
C
C------------------------------------------------------------------------
C
C
C             ...handle special case of isolated atom(s) only.
C                This requires straightforward consecutive indexing.
C
C
         IF (N2CEN.EQ.0) THEN
             DO ATOM = 1,NATOM
                ATORD (ATOM) = ATOM
             END DO
             RETURN
         END IF
C
C
C             ...initialize used atoms index array.
C
C
         USED = 1
         NOTUSED = 0

         DO ATOM = 1,NATOM
            INDEX (ATOM) = NOTUSED
         END DO
C
C
C             ...initialize atom index list.
C
C
         M = 1
         IDX = 1
         ATOM = AT2CEN (1,1)
         ATORD (1) = ATOM
         INDEX (ATOM) = USED
C
C
C             ...check presence of atom ATOM in the set of ordered
C                index pairs and add the new atoms to the atom index
C                list.
C
C
 1000    DO I = 1,N2CEN
            AT1 = AT2CEN (1,I)
            AT2 = AT2CEN (2,I)
            CASE1 = AT1 .EQ. ATOM
            CASE2 = AT2 .EQ. ATOM
            CHECK = CASE1 .OR. CASE2

            IF (CHECK) THEN
                IF (CASE1) THEN
                    ATNEW = AT2
                ELSE
                    ATNEW = AT1
                END IF

                IGNORE = .FALSE.
                DO J = 1,M
                   IGNORE = IGNORE .OR. (ATORD (J).EQ.ATNEW)
                END DO

                IF (.NOT.IGNORE) THEN
                    M = M + 1
                    ATORD (M) = ATNEW
                    INDEX (ATNEW) = USED
                END IF
            END IF
         END DO

         IDX = IDX + 1
         IF (IDX.LE.M) THEN
C
C
C             ...more atoms are present in current atom index list
C                ready for checking. Pick the next in line and go
C                back.
C
C
             ATOM = ATORD (IDX)
             GOTO 1000
         ELSE
C
C
C             ...all atoms in current atom index list have been
C                checked. Check next, if still atoms need to be
C                added to this list corresponding to untreated
C                index pairs. If so, take next available atom and
C                go back to the checking procedure.
C
C
             DO I = 1,N2CEN
                AT1 = AT2CEN (1,I)
                AT2 = AT2CEN (2,I)
                CASE1 = INDEX (AT1) .EQ. NOTUSED
                CASE2 = INDEX (AT2) .EQ. NOTUSED
                CHECK = CASE1 .OR. CASE2

                IF (CHECK) THEN
                    IF (CASE1) THEN
                        ATOM = AT1
                    ELSE
                        ATOM = AT2
                    END IF
                    M = IDX
                    ATORD (IDX) = ATOM
                    INDEX (ATOM) = USED
                    GOTO 1000
                END IF
             END DO
C
C
C             ...we get here if all atoms of all index pairs
C                have been ordered. There might be still atoms
C                which are not bonded to any of the rest and
C                which still must be added to the atom index list
C                in order to complete it.
C
C
             DO ATOM = 1,NATOM
                IF (INDEX (ATOM).EQ.NOTUSED) THEN
                    M = M + 1
                    ATORD (M) = ATOM
                END IF
             END DO
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
