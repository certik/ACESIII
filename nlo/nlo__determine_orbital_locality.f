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
         SUBROUTINE  NLO__DETERMINE_ORBITAL_LOCALITY
     +
     +                    ( NBAS,NATOM,
     +                      ATIDX,
     +                      BASBEG,BASEND,
     +                      C,
     +                      SHALF,
     +
     +                              LOCAL )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__DETERMINE_ORBITAL_LOCALITY
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine analyzes the locality of the orbitals,
C                which are given as the coefficient matrix in AO basis.
C                The quantity that is determined here is the atomic
C                localization, which measures the degree of atomic
C                orbital involvment in each orbital. We can extract
C                the content of the correponding AOs in the coefficent
C                expansion via the following analysis:
C
C                In order to get rid of the AO nonorthogonality
C                problem, we express each orbital in terms of the
C                orthonormal symmetrized AO basis (SAO), which
C                resembles most closely the original AO basis. The
C                new orbital expansion coefficients c' in terms of SAOs
C                are then given by:
C
C                                 c' = S**(1/2) * c
C
C                where c is the coefficient vector for the orbitals in
C                AO basis. Suppose that we want to calculate the atomic
C                content due to atom A in an orbital. Then we form
C                the quantity:
C
C                                 x = c'(T)*c'
C
C                where we restrict the summation over elements of c'
C                only to those corresponding to atom A. This in turn
C                means that we only take those rows of S**(1/2)
C                corresponding to atom A in evaluating c'. The obtained
C                number x we identify as the atomic localization
C                content corresponding to atom A in the orbital.
C
C                The routine evaluates the complete atomic localization
C                map LOCAL (ATOM,NXO) for each atom ATOM and orbital
C                NXO.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    ATIDX (A)    =  will contain atomic index for
C                                    atom A.
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
C                    C            =  NBAS x NBAS orbital coefficient
C                                    matrix in AO basis.
C                    SHALF        =  full NBAS x NBAS square root of
C                                    the overlap matrix in AO basis.
C
C
C                  Output:
C
C                    LOCAL (A,I)  =  atomic localization content of
C                                    atom A for I-th orbital.
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

         LOGICAL     TOTAL

         INTEGER     ATOM
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NXO

         INTEGER     ATIDX  (1:NATOM)
         INTEGER     BASBEG (1:NATOM)
         INTEGER     BASEND (1:NATOM)

         DOUBLE PRECISION  C     (1:NBAS ,1:NBAS)
         DOUBLE PRECISION  LOCAL (1:NATOM,1:NBAS)
         DOUBLE PRECISION  SHALF (1:NBAS ,1:NBAS)
C
C
C------------------------------------------------------------------------
C
C
C             ...form general atomic index vector.
C
C
         DO ATOM = 1,NATOM
            ATIDX (ATOM) = ATOM
         END DO
C
C
C             ...calculate the individual atomic localization contents
C                for each orbital.
C
C
         TOTAL = .FALSE.

         DO NXO = 1,NBAS
            CALL  NLO__EXTRACT_ATOMIC_CONTENT
     +
     +                 ( NBAS,NATOM,
     +                   NATOM,
     +                   ATIDX,
     +                   BASBEG,BASEND,
     +                   TOTAL,
     +                   C (1,NXO),
     +                   SHALF,
     +
     +                           LOCAL (1,NXO) )
     +
     +
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
