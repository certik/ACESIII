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
         LOGICAL FUNCTION  NLO__ACCEPT_ATOMS_BOND_SEARCH
     +
     +                    ( NATOM,
     +                      BONDSIZE,
     +                      NHCEN,NHA,
     +                      ATHIDX,ATHVAL,
     +                      NHYB,
     +                      BOMAT,
     +                      CHOOSED,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      IVEC )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ACCEPT_ATOMS_BOND_SEARCH
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine checks, if a set of atomic indices will
C                be accepted for further bond formation search. There
C                are three resons why the present set of atomic indices
C                might be rejected:
C
C                     1) At least one of the atoms has already a
C                        complete pre-NHO set.
C
C                     2) The set contains isolated (not bonded) atoms,
C                        as judged by the simplified bond order matrix.
C
C                     3) The present particular set of atoms has
C                        already been searched for bond formation
C                        through the bond choosing.
C
C                  Input:
C
C                    NATOM        =  total # of atomic centers
C                    NHCEN        =  current # of atomic hybrid centers
C                                    to be checked for NHCEN centered
C                                    bond construction.
C                    NHA          =  total # of hybrid atoms.
C                    ATHIDX       =  contains the NHCEN atomic indices
C                                    to be tested for bond formation.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    NHYB (A)     =  current total # of pre-NHO's on
C                                    each hybrid atom A.
C                    BOMAT (A,B)  =  simplified bond order matrix
C                                    containing info if atoms A and B
C                                    are considered to be bonded or not.
C                    CHOOSED      =  is true, if at least one NHCEN
C                                    atomic indicex combination has
C                                    been chosen for bond search
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
         LOGICAL     BONDED
         LOGICAL     CHOOSED
         LOGICAL     EQUAL
         LOGICAL     INCRESE

         INTEGER     ATOMI,ATOMJ
         INTEGER     BONDSIZE
         INTEGER     I,J,N
         INTEGER     MXCHOOSE
         INTEGER     NATOM
         INTEGER     NCHOOSE
         INTEGER     NHCEN
         INTEGER     NHA

         INTEGER     ATHIDX  (1:NHCEN)
         INTEGER     ATHVAL  (1:NHA  )
         INTEGER     IVEC    (1:2*NHCEN)
         INTEGER     NHYB    (1:NHA  )

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE)

         DOUBLE PRECISION  ONE

         DOUBLE PRECISION  BOMAT (1:NATOM,1:NATOM)

         DATA  ONE  /1.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...check for complete pre-NHO sets on all atoms
C                involved.
C
C
         ACCEPT = .TRUE.

         DO I = 1,NHCEN
            ATOMI = ATHIDX (I)
            IF (NHYB (ATOMI) .EQ. ATHVAL (ATOMI)) THEN
                ACCEPT = .FALSE.
            END IF
         END DO
C
C
C             ...check for isolated nonbonded atoms.
C
C
         IF (ACCEPT) THEN
             DO J = 1,NHCEN
                ATOMJ = ATHIDX (J)
                BONDED = .FALSE.
                DO I = 1,NHCEN
                   ATOMI = ATHIDX (I)
                   IF (BOMAT (ATOMI,ATOMJ) .EQ. ONE) THEN
                       BONDED = .TRUE.
                   END IF
                END DO
                IF (.NOT.BONDED) THEN
                    ACCEPT = .FALSE.
                END IF
             END DO
         END IF
C
C
C             ...check for already choosen atom combinations.
C
C
         IF (ACCEPT) THEN
             IF (CHOOSED) THEN

                 ABSOLUT = .FALSE.
                 INCRESE = .TRUE.

                 CALL  NLO__SORT_INT_VECTOR_ELEMENTS
     +
     +                      ( NHCEN,NHCEN,
     +                        1,NHCEN,
     +                        ABSOLUT,INCRESE,
     +                        0,
     +                        ATHIDX,
     +
     +                                 IVEC )
     +
     +
                 N = IABS (NCHOOSE)

                 DO I = 1,N

                    CALL  NLO__SORT_INT_VECTOR_ELEMENTS
     +
     +                         ( NHCEN,NHCEN,
     +                           1,NHCEN,
     +                           ABSOLUT,INCRESE,
     +                           0,
     +                           CHOOSE (1,I),
     +
     +                                    IVEC (NHCEN+1))
     +
     +
                    EQUAL = .TRUE.
                    DO J = 1,NHCEN
                       ATOMI = ATHIDX(IVEC (J))
                       ATOMJ = CHOOSE(IVEC (NHCEN+J),I)
                       IF (ATOMI.NE.ATOMJ) THEN
                           EQUAL = .FALSE.
                       END IF
                    END DO

                    IF (EQUAL) THEN
                        NLO__ACCEPT_ATOMS_BOND_SEARCH = .FALSE.
                        RETURN
                    END IF

                 END DO

             END IF
         END IF

         NLO__ACCEPT_ATOMS_BOND_SEARCH = ACCEPT
C
C
C             ...ready!
C
C
         RETURN
         END
