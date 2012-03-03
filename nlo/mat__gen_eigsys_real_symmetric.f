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
         SUBROUTINE  MAT__GEN_EIGSYS_REAL_SYMMETRIC
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWS,DDCOLS,
     +                      DDVECD,
     +                      N,
     +                      REVERS,
     +
     +                              D,
     +                              S,
     +                              X )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__GEN_EIGSYS_REAL_SYMMETRIC
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : Computation of all eigenvalues and corresponding
C                eigenvectors of a real symmetric generalized
C                eigenvalue problem:
C
C
C                             Ax = cSx
C
C
C                where A is symmetric and S is symmetric and positive
C                definite.
C
C
C                Procedure:
C                ----------
C
C                The matrix S is decomposed into: S = B * B (T), where
C                B is a lower triangular matrix and (T) denotes
C                transposition. Then the above generalized eigenvalue
C                problem can be rewritten as:
C
C
C                   {B(inv) A B(inv,T)} {B(T) x}   =   c {B(T) x}
C
C
C                where (inv) denotes inversion. Hence we obtained
C                an ordinary eigenvalue problem, which can be solved
C                by a conventional technique.
C
C
C                  Input:
C
C                   DDROWz = row dimensions for matrices X and S
C                            (z = X,S)
C                   DDCOLz = column dimensions for matrices X and S
C                            (z = X,S)
C                   DDVECD = dimension for vector D containing the
C                            eigenvalues
C                        N = order of the matrix X to be diagonalized
C                        X = matrix to be diagonalized
C                        S = overlap matrix (must be positive definite)
C                   REVERS = if false, eigenvalues and associated
C                            eigenvectors are in ascending sequence,
C                            otherwise this ordering is reversed.
C
C                  Output:
C
C                        D = computed eigenvalues, in specified order
C                        X = matrix of corresponding eigenvectors,
C                            columnwise.
C
C                Only the lower triangle of X and S need to be
C                supplied. However, both X and S will be destroyed
C                during the diagonalization process!!!
C
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         LOGICAL   REVERS

         INTEGER   DDROWX,DDCOLX,DDROWS,DDCOLS,DDVECD
         INTEGER   I,J,L,N

         DOUBLE PRECISION   DELTA
         DOUBLE PRECISION   ONE
         DOUBLE PRECISION   ROOT
         DOUBLE PRECISION   SII
         DOUBLE PRECISION   SQRARG
         DOUBLE PRECISION   SUM
         DOUBLE PRECISION   ZERO

         DOUBLE PRECISION   D (1:DDVECD)

         DOUBLE PRECISION   S (1:DDROWS,1:DDCOLS)
         DOUBLE PRECISION   X (1:DDROWX,1:DDCOLX)

         DATA   DELTA    /1.D-17/
         DATA   ZERO     /0.D0 /
         DATA   ONE      /1.D0 /
C
C
C------------------------------------------------------------------------
C
C
C             ...check passed dimensions.
C
C
         IF ( N.GT.DDROWX .OR. N.GT.DDCOLX )  THEN
              WRITE (1,*) ' Dimension of matrix X too small: '
              WRITE (1,*) ' mat__gen_eigsys_real_symmetric '
              WRITE (1,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
              WRITE (*,*) ' Dimension of matrix X too small: '
              WRITE (*,*) ' mat__gen_eigsys_real_symmetric '
              WRITE (*,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
              STOP
         END IF

         IF ( N.GT.DDROWS .OR. N.GT.DDCOLS )  THEN
              WRITE (1,*) ' Dimension of overlap matrix S too small: '
              WRITE (1,*) ' mat__gen_eigsys_real_symmetric '
              WRITE (1,*) ' DDROWS,DDCOLS,N = ',DDROWS,DDCOLS,N
              WRITE (*,*) ' Dimension of overlap matrix S too small: '
              WRITE (*,*) ' mat__gen_eigsys_real_symmetric '
              WRITE (*,*) ' DDROWS,DDCOLS,N = ',DDROWS,DDCOLS,N
              STOP
         END IF

         IF ( N.GT.DDVECD )  THEN
              WRITE (1,*) ' Dimension of vector D too small: '
              WRITE (1,*) ' mat__gen_eigsys_real_symmetric '
              WRITE (1,*) ' DDVECD,N = ',DDVECD,N
              WRITE (*,*) ' Dimension of vector D too small: '
              WRITE (*,*) ' mat__gen_eigsys_real_symmetric '
              WRITE (*,*) ' DDVECD,N = ',DDVECD,N
              STOP
         END IF
C
C
C             ...determine B from S via a Cholesky decomposition
C                and place B in upper part of S array.
C
C
          IF ( S (1,1).LE.DELTA ) THEN
               WRITE (*,*) ' S matrix not (enough) +ve definite! '
               WRITE (*,*) ' mat__gen_eigsys_real_symmetric '
               WRITE (*,*) ' Square root argument = ',S (1,1)
               WRITE (1,*) ' S matrix not (enough) +ve definite! '
               WRITE (1,*) ' mat__gen_eigsys_real_symmetric '
               WRITE (1,*) ' Square root argument = ',S (1,1)
               STOP
          END IF

          ROOT  =  DSQRT ( S (1,1) )
          S (1,1)  =  ROOT

          DO  10  I = 2,N
              S (1,I)  =  S (I,1) / ROOT
   10     CONTINUE

          DO  20  J = 2,N

              SQRARG  =  S (J,J)
              DO  30  L = 1,J-1
                  SQRARG  =  SQRARG  -  S (L,J) ** 2
   30         CONTINUE

              IF ( SQRARG.LE.DELTA ) THEN
                   WRITE (*,*) ' S matrix not (enough) +ve definite! '
                   WRITE (*,*) ' mat__gen_eigsys_real_symmetric '
                   WRITE (*,*) ' Square root argument = ',SQRARG
                   WRITE (1,*) ' S matrix not (enough) +ve definite! '
                   WRITE (1,*) ' mat__gen_eigsys_real_symmetric '
                   WRITE (1,*) ' Square root argument = ',SQRARG
                   STOP
              END IF

              ROOT  =  DSQRT (SQRARG)
              S (J,J)  =  ROOT

              DO  40  I = J+1,N
                  SUM  =  S (I,J)
                  DO  50  L = 1,J-1
                      SUM  =  SUM  -  S (L,I) * S (L,J)
   50             CONTINUE
                  S (J,I)  =  SUM / ROOT
   40         CONTINUE

   20     CONTINUE
C
C
C
C             ...calculate next the product:  C = B(inv) A B(inv,T).
C                First observe that, since the result is symmetric,
C                we need to calculate only the upper triangle of C.
C                From this follows, since B(inv) is lower triangular,
C                that one only needs to calculate the lower triangle
C                of the product: B(inv) A. Hence, the following two
C                steps are performed:
C
C                  1) Form lower triangle of:  B(inv) A  and store
C                     the result in the upper part of X
C                
C                  2) Premultiply the result of 1) by B(inv) by
C                     overwriting the upper part of X
C                
C                After 1) and 2) the upper triangle of C calculated
C                sits in the upper triangle of X.
C
C                Note further, that actual calculation of B(inv) is
C                not needed. This matrix is 'calculated' in an
C                implicit way from B, using its nice lower triangular
C                property.
C
C
C
          DO  100  I = 1,N
              SII  =  S (I,I)
              DO  200  J = I,N
                  SUM  =  X (J,I)
                  DO  300  L = 1,I-1
                      SUM  =  SUM  -  S (L,I) * X (L,J)
  300             CONTINUE
                  X (I,J)  =  SUM / SII
  200         CONTINUE
  100     CONTINUE

          DO  400  J = 1,N
          DO  400  I = J,N
              SUM  =  X (J,I)
              DO  500  L = J,I-1
                  SUM  =  SUM  -  X (J,L) * S (L,I)
  500         CONTINUE
              DO  600  L = 1,J-1
                  SUM  =  SUM  -  X (L,J) * S (L,I)
  600         CONTINUE
              X (J,I)  =  SUM / S (I,I)
  400     CONTINUE
C
C
C             ...move upper triangle of X to lower triangle of X
C                to be able to call the standard diagonalization
C                routine.
C
C
          DO  700  I = 1,N
          DO  700  J = I,N
              X (J,I) = X (I,J)
  700     CONTINUE

          CALL    MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                 ( DDROWX,DDCOLX,DDVECD,
     +                   N,
     +                   REVERS,
     +
     +                        D,
     +                        X )
     +
     +
C
C
C             ...transform back the eigenvectors of the standard case
C                to the eigenvectors of the general case by making the
C                product:  B(inv,T) X  and overwriting X.
C
C
          DO  800  J = 1,N
          DO  800  I = N,1,-1
              SUM  =  X (I,J)
              DO  900  L = I+1,N
                  SUM  =  SUM  -  S (I,L) * X (L,J)
  900         CONTINUE
              X (I,J)  =  SUM / S (I,I)
  800     CONTINUE
C
C
C             ...finished!
C
C
          RETURN
          END
