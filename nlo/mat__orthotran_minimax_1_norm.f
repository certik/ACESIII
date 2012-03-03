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
         SUBROUTINE  MAT__ORTHOTRAN_MINIMAX_1_NORM
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWT,DDCOLT,
     +                      DDW1,DDW2,
     +                      ROW,COL,
     +                      NORMLZE,
     +                      MAXIMZE,
     +                      X,
     +                      W1,W2,
     +
     +                              NITER,
     +                              T )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__ORTHOTRAN_MINIMAX_1_NORM
C  MODULE      : Matrix
C  MODULE-ID   : MAT
C  DESCRIPTION : This routine finds the orthogonal rotation matrix T
C                that will lead to a minimization or maximization of
C                the 1-norm of a set of vectors X:
C
C                         X' = X * T     ||X'||  = minimum or maximum
C                                              1
C                Theory:
C                ------
C
C                Consider first a pair of vectors X1 and X2. Then we
C                need to find the orthogonal 2x2 rotation matrix given
C                by the two rotation angles C = cos (theta) and
C                S = sin (theta) with condition C*C + S*S = 1:
C
C                                      C   S
C
C                                     -S   C
C
C                         (X1   X2)  (X1'  X2')
C
C                Thus we need to find either C or S that minimize or
C                maximize the expression for the overall 1-norm:
C
C                       ||C*X1 - S*X2||  + ||S*X1 + C*X2||
C                                      1                  1
C
C                Maximization of this expression will lead (after
C                differentiation) to a discontinous function in X1
C                and X2 components possessing 1st derivative jumps.
C                Hence we reformulate the problem and, assuming the
C                X1 and X2 are normalized (i.e. all their components
C                =< 1), we can state the problem as maximizing or
C                minimizing the sum of all 4th powers of the resulting
C                rotated vector components in X1' and X2':
C
C                     sum    (C * X1 (k) - S * X2 (k)) ** 4
C                      k   + (S * X1 (k) + C * X2 (k)) ** 4
C
C                When differentiating this sum with respect to either
C                C or S and making the variable substitution Z = S / C
C                we obtain after some manipulation a 4th order
C                polynomial equation in Z:
C
C                   Z**4 - (B/A) Z**3 - 6 (Z**2) + (B/A)*Z + 1 = 0
C
C                where the two values of A and B are given by:
C
C                 A  =  sum  X2(k) * X1(k)**3 - X1(k) * X2(k)**3
C                        k
C
C                 B  =  sum  X1(k)**4 + X2(k)**4 - 6 * (X1(k)*X2(k))**2
C                        k
C
C                which can be factored into two quadratic equations:
C
C                    (Z**2 + P * Z - 1) (Z**2 + Q * Z - 1) = 0
C
C                Multiplying the factors out and equating coefficients
C                with the 4-th order polynomial equation in Z we arrive
C                at the two equations in P and Q:
C
C                              P + Q = - (B/A)
C
C                                 PQ = - 4
C
C                leading to quadratic equations with the P,Q pair
C                solution:
C
C                             - (B/A) + SQRT ((B/A)**2 + 16)
C                         P = ------------------------------
C                                           2
C
C                             - (B/A) - SQRT ((B/A)**2 + 16)
C                         Q = ------------------------------
C                                           2
C
C                A total of 4 solutions of Z are thus possible, of
C                which 2 will characterize the maximum and minimum
C                rotation angle. The first two solutions come from
C                the P-quadratic equation, giving:
C
C
C
C                Note that the absolute value of A determines how
C                large the rotation angle will be. If |A| is very
C                close to zero, then the two vectors X1 and X2 are
C                already very close to their minimum or maximum 1-norm
C                form.
C
C                Procedure:
C                ---------
C
C                If more than two vectors are present in X, i.e. if
C                we have COL > 2, we adopt the following procedure:
C
C                   1) Determine the pair of vectors which give the
C                      largest |A| value.
C
C                   2) If largest |A| value < threshold, we are done
C                      => exit.
C
C                   3) Find 2x2 rotation matrix elements C and S
C                      that minimizes the sum of the 4-th powers
C                      according to the theory outlined above.
C
C                   4) Accumulate C and S to the final COL x COL
C                      orthogonal rotation matrix T.
C
C                   5) Return to 1)
C
C
C
C                  Input:
C
C                    DDROWX,DDCOLX  =  declared dimensions of vectors X.
C                    DDROWT,DDCOLT  =  declared dimensions of rotation
C                                      matrix T.
C                    DDW1,DDW2      =  declared dimensions of working
C                                      arrays W1 and W2.
C                    ROW,COL        =  actual dimensions of vectors X.
C                    NORMLZE        =  if true, the initial vectors X
C                                      will be normalized.
C                    MAXIMZE        =  if true, the 1-norm will be
C                                      maximized, otherwise the 1-norm
C                                      will be minimized.
C                    X              =  initial ROW x COL vectors X.
C                    W1,W2          =  working arrays.
C
C
C                  Output:
C
C                    NITER          =  # of iterations performed.
C                    T              =  orthogonal rotation matrix.
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
C         IMPLICIT    DOUBLE PRECISION (A-H,O-Z)

         LOGICAL     MAXIMZE
         LOGICAL     NORMLZE
         LOGICAL     ROTATE

         INTEGER     DDROWT,DDCOLT
         INTEGER     DDROWX,DDCOLX
         INTEGER     DDW1,DDW2
         INTEGER     I,J,K,L,N
         INTEGER     NITER,MXITER
         INTEGER     ROW,COL

         DOUBLE PRECISION  A,B,C,D,R
         DOUBLE PRECISION  AMAX
         DOUBLE PRECISION  C1,S1,C2,S2
         DOUBLE PRECISION  ROOT1,ROOT2
         DOUBLE PRECISION  SMALL,VSMALL,CONVGED
         DOUBLE PRECISION  SUM1,SUM2
         DOUBLE PRECISION  XKK,XKL,XLL
         DOUBLE PRECISION  ZERO,HALF,ONE,TWO,FOUR,SIX,SIXTEEN

         DOUBLE PRECISION  W1 (1:DDW1)
         DOUBLE PRECISION  W2 (1:DDW2)

         DOUBLE PRECISION  T (1:DDROWT,1:DDCOLT)
         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)

         PARAMETER  (ZERO    = 0.D0)
         PARAMETER  (ONE     = 1.D0)
         PARAMETER  (TWO     = 2.D0)
         PARAMETER  (FOUR    = 4.D0)
         PARAMETER  (SIX     = 6.D0)
         PARAMETER  (SIXTEEN = 16.D0)
         PARAMETER  (HALF    = 0.5D0)
         PARAMETER  (MXITER  = 1000)
         PARAMETER  (SMALL   = 1.D-6)
         PARAMETER  (VSMALL  = 1.D-9)
         PARAMETER  (CONVGED = 1.D-14)
C
C
C------------------------------------------------------------------------
C
C
C             ...return immediately, if COL = 1.
C
C
         IF (COL.EQ.1) THEN
             NITER = 1
             T (1,1) = ONE
             RETURN
         END IF
C
C
C             ...check dimensions for:
C
C                      i) vectors X
C                     ii) final rotation matrix T
C                    iii) working vectors W1 and W2
C
C
         IF  (ROW.GT.DDROWX .OR. COL.GT.DDCOLX) THEN
              WRITE (*,*) ' Dimensions of vectors X too small! '
              WRITE (*,*) ' mat__orthotran_minimax_1_norm '
              WRITE (*,*) ' DDROWX,DDCOLX,ROW,COL = ',
     +                      DDROWX,DDCOLX,ROW,COL
              WRITE (1,*) ' Dimensions of vectors X too small! '
              WRITE (1,*) ' mat__orthotran_minimax_1_norm '
              WRITE (1,*) ' DDROWX,DDCOLX,ROW,COL = ',
     +                      DDROWX,DDCOLX,ROW,COL
              STOP
         END IF

         IF  (COL.GT.DDROWT .OR. COL.GT.DDCOLT) THEN
              WRITE (*,*) ' Dimensions of rot matrix T too small! '
              WRITE (*,*) ' mat__orthotran_minimax_1_norm '
              WRITE (*,*) ' DDROWT,DDCOLT,COL = ',DDROWT,DDCOLT,COL
              WRITE (1,*) ' Dimensions of rot matrix T too small! '
              WRITE (1,*) ' mat__orthotran_minimax_1_norm '
              WRITE (1,*) ' DDROWT,DDCOLT,COL = ',DDROWT,DDCOLT,COL
              STOP
         END IF

         IF  (ROW.GT.DDW1 .OR. ROW.GT.DDW2) THEN
              WRITE (*,*) ' Dimensions of vectors W1 or W2 too small! '
              WRITE (*,*) ' mat__orthotran_minimax_1_norm '
              WRITE (*,*) ' DDW1,DDW2,ROW = ',DDW1,DDW2,ROW
              WRITE (1,*) ' Dimensions of vectors W1 or W2 too small! '
              WRITE (1,*) ' mat__orthotran_minimax_1_norm '
              WRITE (1,*) ' DDW1,DDW2,ROW = ',DDW1,DDW2,ROW
              STOP
         END IF
C
C
C             ...normalize the vectors in X (if wanted). If the maximum
C                norm difference is greater than 10-4 percent of the
C                minimum norm, the routine issues a warning! 
C
C
         IF (NORMLZE) THEN

             B = ZERO
             DO N = 1,ROW
                B = B + X (N,1) * X (N,1)
             END DO
             B = ONE / DSQRT (B)
             DO N = 1,ROW
                X (N,1) = B * X (N,1)
             END DO

             A = ZERO
             C = ZERO
             DO K = 2,COL
                D = ZERO
                DO N = 1,ROW
                   D = D + X (N,K) * X (N,K)
                END DO
                D = ONE / DSQRT (D)
                DO N = 1,ROW
                   X (N,K) = D * X (N,K)
                END DO
                A = MIN (A,D-B)
                C = MAX (C,D-B)
             END DO
             D = B + A
             B = C - A
             A = B / D

             IF (A.GT.1.D-6) THEN
                 WRITE (*,*) ' WARNING! Norms too different! '
                 WRITE (*,*) ' mat__orthotran_minimax_1_norm '
                 WRITE (*,*) ' Minimum norm            = ',D
                 WRITE (*,*) ' Maximum norm difference = ',B
                 WRITE (*,*) ' Proceeding ... '
                 WRITE (1,*) ' WARNING! Norms too different! '
                 WRITE (1,*) ' mat__orthotran_minimax_1_norm '
                 WRITE (1,*) ' Minimum norm            = ',D
                 WRITE (1,*) ' Maximum norm difference = ',B
                 WRITE (1,*) ' Proceeding ... '
             END IF

         END IF
C
C
C             ...initialize rotation matrix T.
C
C
         DO J = 1,COL
            DO I = 1,COL
               T (I,J) = ZERO
            END DO
            T (J,J) = ONE
         END DO
C
C
C             ...prepare for small perturbed rotation on pair of
C                vectors possessing already converged value of |A|.
C
C
         C1 = DSQRT (ONE - SMALL)
         S1 = DSQRT (SMALL)

         DO K = 1,COL
         DO L = K+1,COL
            A = ZERO
            B = ZERO
            C = ZERO
            D = ZERO
            DO N = 1,ROW
               XKK = X (N,K) * X (N,K)
               XKL = X (N,K) * X (N,L)
               XLL = X (N,L) * X (N,L)
               A = A + XKK * XKL
               B = B + XKK * XKK + XLL * XLL
               C = C + XKL * XLL
               D = D + XKL * XKL
            END DO
            A = A - C
            B = B - SIX * D

            WRITE (*,*) ' A , B (seed) = ',A,B
            ROTATE = DABS (A).LT.CONVGED .AND.
     +               DABS (B).GT.VSMALL

            IF (ROTATE) THEN
                DO N = 1,ROW
                   A = C1 * X (N,K) - S1 * X (N,L)
                   C = S1 * X (N,K) + C1 * X (N,L)
                   X (N,K) = A
                   X (N,L) = C
                END DO
                DO N = 1,COL
                   A = C1 * T (N,K) - S1 * T (N,L)
                   C = S1 * T (N,K) + C1 * T (N,L)
                   T (N,K) = A
                   T (N,L) = C
                END DO
            END IF
         END DO
         END DO
C
C
C             ...perform iterates.
C
C
         DO NITER = 1,MXITER

            AMAX = CONVGED
C
C
C             ...find I,J pair of X vectors possessing largest |A|
C                value.
C
C
            DO K = 1,COL
            DO L = K+1,COL
               A = ZERO
               B = ZERO
               C = ZERO
               D = ZERO
               DO N = 1,ROW
                  XKK = X (N,K) * X (N,K)
                  XKL = X (N,K) * X (N,L)
                  XLL = X (N,L) * X (N,L)
                  A = A + XKK * XKL
                  B = B + XKK * XKK + XLL * XLL
                  C = C + XKL * XLL
                  D = D + XKL * XKL
               END DO
               A = A - C
               B = B - SIX * D
               WRITE (*,*) ' A , B (iter) = ',A,B
               C = DABS (A)
               IF (C.GT.AMAX) THEN
                   I = K
                   J = L
                   R = B / A
                   AMAX = C
               END IF
            END DO
            END DO

            WRITE (*,*) ' AMAX = ',AMAX
C
C
C             ...check, if convergence was achieved (largest |A|
C                below convergence threshold) and exit if the case.
C
C
            IF (AMAX.EQ.CONVGED) THEN
                WRITE (*,*) ' Converged after ',NITER,' iterations! '
                CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                     ( 6,
     +                       ' Final set of vectors ',
     +                       DDROWX,DDCOLX,
     +                       ROW,COL,
     +                       X )
     +
c                CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
c     +
c     +                     ( 6,
c     +                       ' Rotation matrix ',
c     +                       DDROWT,DDCOLT,
c     +                       COL,COL,
c     +                       T )
c     +
                RETURN
            END IF

c            WRITE (*,*) ' Processing columns ',I,J
C
C
C             ...start rotation of I,J vector pair. Calculate the
C                2 minimax roots.
C
C
            IF (R.GE.ZERO) THEN
                ROOT2 = - HALF * (R + DSQRT (R*R + SIXTEEN))
                ROOT1 = - ONE / ROOT2
            ELSE
                ROOT1 = - HALF * (R - DSQRT (R*R + SIXTEEN))
                ROOT2 = - ONE / ROOT1
            END IF

c            WRITE (*,*) ' Root 1 , Root 2 = ',ROOT1,ROOT2
C
C
C             ...calculate both sets of C1,S1 and C2,S2 values for
C                both roots ROOT1 and ROOT2, avoiding loss of accuracy
C                as much as possible for small angles.
C
C
            A = ONE / (ROOT1 * ROOT1 + FOUR)
            B = ONE / (ROOT2 * ROOT2 + FOUR)

            IF (A.LT.SMALL) THEN
                A = A + A * A
            ELSE
                A = (ONE - DSQRT (ONE - FOUR * A)) * HALF
            END IF

            IF (B.LT.SMALL) THEN
                B = B + B * B
            ELSE
                B = (ONE - DSQRT (ONE - FOUR * B)) * HALF
            END IF

            S1 = DSQRT (A)
            C1 = DSQRT (ONE - A)
            C2 = DSQRT (B)
            S2 = DSQRT (ONE - B)
C
C
C             ...calculate both 4th power summations for both roots.
C                Near the middle to the end of convergence, the largest
C                of the absolute root values is the one leading to the
C                right direction of minimization or maximization. Hence
C                the rotated I,J vector pair are overwritten by the
C                rotated results which correspond to the largest
C                absolute root value. Only in case (during the beginning
C                few iterations) the smallest absolute root value is
C                the one leading to the right direction we copy the
C                rotated vector pair results stored in W1 and W2 for
C                the smallest root back to the original vector pair
C                memory.
C
C
            SUM1 = ZERO
            SUM2 = ZERO

            IF (DABS (ROOT1) .GT. DABS (ROOT2)) THEN

                DO N = 1,ROW
                   A = C1 * X (N,I) - S1 * X (N,J)
                   B = S1 * X (N,I) + C1 * X (N,J)
                   C = C2 * X (N,I) - S2 * X (N,J)
                   D = S2 * X (N,I) + C2 * X (N,J)
                   X (N,I) = A
                   X (N,J) = B
                   W1 (N)  = C
                   W2 (N)  = D
                   A = A * A
                   B = B * B
                   C = C * C
                   D = D * D
                   SUM1 = SUM1 + (A * A) + (B * B)
                   SUM2 = SUM2 + (C * C) + (D * D)
                END DO

                IF (MAXIMZE) THEN
                    IF (SUM2.LT.SUM1) THEN
                        C1 = C2
                        S1 = S2
                        DO N = 1,ROW
                           X (N,I) = W1 (N)
                           X (N,J) = W2 (N)
                        END DO
                    END IF
                ELSE
                    IF (SUM2.GT.SUM1) THEN
                        C1 = C2
                        S1 = S2
                        DO N = 1,ROW
                           X (N,I) = W1 (N)
                           X (N,J) = W2 (N)
                        END DO
                    END IF
                END IF

            ELSE

                DO N = 1,ROW
                   A = C1 * X (N,I) - S1 * X (N,J)
                   B = S1 * X (N,I) + C1 * X (N,J)
                   C = C2 * X (N,I) - S2 * X (N,J)
                   D = S2 * X (N,I) + C2 * X (N,J)
                   X (N,I) = C
                   X (N,J) = D
                   W1 (N)  = A
                   W2 (N)  = B
                   A = A * A
                   B = B * B
                   C = C * C
                   D = D * D
                   SUM1 = SUM1 + (A * A) + (B * B)
                   SUM2 = SUM2 + (C * C) + (D * D)
                END DO

                IF (MAXIMZE) THEN
                    IF (SUM1.LT.SUM2) THEN
                        DO N = 1,ROW
                           X (N,I) = W1 (N)
                           X (N,J) = W2 (N)
                        END DO
                    ELSE
                        C1 = C2
                        S1 = S2
                    END IF
                ELSE
                    IF (SUM1.GT.SUM2) THEN
                        DO N = 1,ROW
                           X (N,I) = W1 (N)
                           X (N,J) = W2 (N)
                        END DO
                    ELSE
                        C1 = C2
                        S1 = S2
                    END IF
                END IF

            END IF

c            WRITE (*,*) ' Sum 1 , Sum 2 = ',SUM1,SUM2
c            CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
c     +
c     +                 ( 6,
c     +                   ' Intermediate set of vectors ',
c     +                   DDROWX,DDCOLX,
c     +                   ROW,COL,
c     +                   X )
c     +
C
C
C             ...update rotation matrix T with proper rotation elements
C                C1 and S1.
C
C
            DO N = 1,COL
               A = C1 * T (N,I) - S1 * T (N,J)
               C = S1 * T (N,I) + C1 * T (N,J)
               T (N,I) = A
               T (N,J) = C
            END DO
C
C
C             ...next iteration.
C
C
         END DO
C
C
C             ...Maximum # of iterations exceeded. Print info.
C
C
         WRITE (*,*) ' Minimax of 1-norm did not converge! '
         WRITE (*,*) ' Maximum # of iterations performed = ',MXITER
         WRITE (*,*) ' mat__orthotran_minimax_1_norm '
         WRITE (1,*) ' Minimax of 1-norm did not converge! '
         WRITE (1,*) ' Maximum # of iterations performed = ',MXITER
         WRITE (1,*) ' mat__orthotran_minimax_1_norm '
C
C
C             ...ready!
C
C
         RETURN
         END
