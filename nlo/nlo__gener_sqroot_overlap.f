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
         SUBROUTINE  NLO__GENER_SQROOT_OVERLAP
     +
     +                    ( NBAS,
     +                      XVEC,XMAT,
     +
     +                              S )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_SQROOT_OVERLAP
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the S**(1/2) matrix and
C                overwrites the original overlap matrix S.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    S            =  full NBAS x NBAS overlap matrix.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C
C
C                  Output:
C
C                    S            =  full NBAS x NBAS S**(1/2) matrix.
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

         INTEGER     I,J,N
         INTEGER     NBAS

         DOUBLE PRECISION  LDEPEND
         DOUBLE PRECISION  X,Y

         DOUBLE PRECISION  XVEC (1:NBAS)

         DOUBLE PRECISION  S    (1:NBAS,1:NBAS)
         DOUBLE PRECISION  XMAT (1:NBAS,1:NBAS)

         DATA  REVERS   /.TRUE./
         DATA  LDEPEND  /1.D-12/
C
C
C------------------------------------------------------------------------
C
C
C             ...generate the S**(1/2) matrix and overwrite original
C                matrix S.
C
C
         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                S,
     +
     +                        XMAT )
     +
     +
         CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +              ( NBAS,NBAS,NBAS,
     +                NBAS,
     +                REVERS,
     +
     +                        XVEC,
     +                        XMAT )
     +
     +
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +
     +                        S )
     +
     +
         DO 100 N = 1,NBAS

            X = XVEC (N)
            IF (X.LT.LDEPEND) THEN
                WRITE (*,*) ' Overlap matrix is linear dependent! '
                WRITE (*,*) ' nlo__gener_sqroot_overlap '
                WRITE (*,*) ' X,LDEPEND = ',X,LDEPEND
                WRITE (1,*) ' Overlap matrix is linear dependent! '
                WRITE (1,*) ' nlo__gener_sqroot_overlap '
                WRITE (1,*) ' X,LDEPEND = ',X,LDEPEND
                STOP
            END IF

            X = DSQRT (X)
            DO 110 J = 1,NBAS
               Y = X * XMAT (J,N)
               DO 120 I = 1,NBAS
                  S (I,J) = S (I,J) + Y * XMAT (I,N)
  120          CONTINUE
  110       CONTINUE

  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
