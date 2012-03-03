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
         SUBROUTINE  NLO__PERTURB_OCCUPATION_MATRIX
     +
     +                    ( DDROWP,DDCOLP,
     +                      NBAS,NATOM,
     +                      ATNBAS,
     +                      EFRAC,
     +                      USED,TEMP,
     +
     +                            P )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PERTURB_OCCUPATION_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine introduces a small perturbation X on
C                the diagonals of the occupation matrix P:
C
C
C                     *                          *
C
C                     *   *                      *   *-X
C
C                     *   *   *         ----->   *   *   *
C
C                     *   *   *   *              *   *   *   *-X
C
C                     *   *   *   *   *          *   *   *   *   *
C
C
C                The positions of the perturbation is determined as
C                follows:
C
C                    Loop over all atomic subblocks of P
C
C                      For each atomic subblock:
C
C                         i) Sum up all diagonals -> sum = S
C                        ii) Determine largest diagonal element -> D
C                       iii) Determine number of diagonal elements
C                            equal to D -> n
C                        iv) Set X = EFRAC * S / n
C                         v) Subtract X from all diagonals equal to D
C
C                This procedure ensures a controled removal of a
C                specific electronic fraction per atomic subblock and
C                at the same time retains any symmetry present in P.
C
C
C                  Input:
C
C                    DDROWP,DDCOLP  =  declared dimensions of matrix P
C                    NBAS           =  dimension of occupation matrix P
C                    NATOM          =  number of atomic center for which
C                                      P is defined
C                    ATNBAS (A)     =  dimension of atomic subblock in
C                                      P corresponding to atom A.
C                    EFRAC          =  electronic fraction to be removed
C                                      from each atomic subblock
C                    USED (I)       =  indicator for which diagonal
C                                      elements of of P have been
C                                      perturbed.
C                    TEMP           =  temporary scratch array holding
C                                      diagonal elements of P
C                    P              =  original occupation matrix
C
C
C                  Output:
C
C                    P              =  perturbed occupation matrix
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

         INTEGER     ATOM
         INTEGER     DDROWP,DDCOLP
         INTEGER     I,J,M
         INTEGER     NATOM
         INTEGER     NBAS,NBASA
         INTEGER     NUSED
         INTEGER     OFFD

         INTEGER     ATNBAS  (1:NATOM)
         INTEGER     USED    (1:NBAS )

         DOUBLE PRECISION  D,DMAX,DTOL
         DOUBLE PRECISION  EFRAC
         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  X
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  TEMP  (1:NBAS)

         DOUBLE PRECISION  P     (1:DDROWP,1:DDCOLP)

         DATA  DTOL  /1.D-10/
         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of P matrix supplied.
C
C
         IF (NBAS.GT.DDROWP .OR. NBAS.GT.DDCOLP) THEN
             WRITE (*,*) ' Dimensions of matrix P too small! '
             WRITE (*,*) ' nlo__perturb_occupation_matrix '
             WRITE (*,*) ' DDROWP,DDCOLP,NBAS = ',DDROWP,DDCOLP,NBAS
             WRITE (1,*) ' Dimensions of matrix P too small! '
             WRITE (1,*) ' nlo__perturb_occupation_matrix '
             WRITE (1,*) ' DDROWP,DDCOLP,NBAS = ',DDROWP,DDCOLP,NBAS
             STOP
         END IF
C
C
C             ...copy all diagonals of P into scratch array.
C
C
         DO 100 I = 1,NBAS
            TEMP (I) = P (I,I)
  100    CONTINUE
C
C
C             ...outer loop over all atomic subblocks.
C
C
         OFFD = 0
         NUSED = 0

         DO 200 ATOM = 1,NATOM

            NBASA = ATNBAS (ATOM)

            SUM = ZERO
            DMAX = ZERO

            DO 210 I = 1,NBASA
               D = TEMP (OFFD+I)
               DMAX = DMAX1 (D,DMAX)
               SUM = SUM + D
  210       CONTINUE

            M = NUSED
            DO 220 I = 1,NBASA
               J = OFFD + I
               D = TEMP (J)
               IF (D.LT.(DMAX+DTOL) .AND. D.GT.(DMAX-DTOL)) THEN
                   NUSED = NUSED + 1
                   USED (NUSED) = J
               END IF
  220       CONTINUE

            X = EFRAC * SUM / DFLOAT (NUSED-M)

            DO 230 I = M+1,NUSED
               J = USED (I)
               TEMP (J) = TEMP (J) - X
  230       CONTINUE

            OFFD = OFFD + NBASA

  200    CONTINUE
C
C
C             ...put perturbed elements of scratch array back into
C                diagonals of P.
C
C
         DO 300 I = 1,NUSED
            J = USED (I)
            P (J,J) = TEMP (J)
  300    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
