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
         SUBROUTINE  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                    ( DDROW,DDCOL,DDVEC,
     +                      N,
     +                      REVERS,
     +
     +                              D,
     +                              X )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__DIAGONALIZE_REAL_SYMMETRIC
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : Computation of all eigenvalues and corresponding
C                eigenvectors of a real symmetric matrix by the method
C                of QL transformations. The code was derived from the
C                routines TRED2 and TQL2 in the EISPACK collection.
C
C                EPS is the relative machine precision of floating
C                point arithmetic, i.e. the minimum of all X such that
C                1.+X is greater than 1. on the computer. If for
C                example a computer has a 48-bit mantissa, then we
C                have EPS = 2**(-47).
C
C                TOL is the smallest positive number for which the sum
C                of the squares of the off-diagonal elements in a row
C                of the (partially transformed) lower triangular
C                matrix is still regarded as nonzero. In the original
C                (slightly different) code this was set to TOL=SPN/EPS,
C                where SPN is the smallest positive number representable
C                within the computer.
C
C                  Input:
C
C                    DDROW = row dimension for matrix X to be
C                            diagonalized
C                    DDCOL = column dimension for matrix X to be
C                            diagonalized
C                    DDVEC = dimension for vector D containing the
C                            eigenvalues
C                        N = order of the matrix X to be diagonalized
C                        X = matrix to be diagonalized
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
C                Note, that the original X is destroyed during the
C                process!!!
C
C
C  AUTHOR      : from EISPACK
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         LOGICAL   REVERS

         INTEGER   DDROW,DDCOL,DDVEC
         INTEGER   I,J,K,L,N
         INTEGER   J1
         INTEGER   LDN
         INTEGER   NI

         PARAMETER    ( LDN = 4000 )

         DOUBLE PRECISION   B,C,F,G,H,P,R,S
         DOUBLE PRECISION   EPS,TOL

         DOUBLE PRECISION   D (1:DDVEC)
         DOUBLE PRECISION   E (1:LDN)

         DOUBLE PRECISION   X (1:DDROW,1:DDCOL)

         DATA   EPS,TOL    /1.D-14,1.D-40/
C
C
C------------------------------------------------------------------------
C
C
C             ...check passed dimension of matrix X and vector D.
C
C
         IF ( N.GT.DDROW .OR. N.GT.DDCOL )  THEN
              WRITE (1,*) ' Dimension of matrix X too small: '
              WRITE (1,*) ' mat__diagonalize_real_symmetric '
              WRITE (1,*) ' DDROW,DDCOL,N = ',DDROW,DDCOL,N
              WRITE (*,*) ' Dimension of matrix X too small: '
              WRITE (*,*) ' mat__diagonalize_real_symmetric '
              WRITE (*,*) ' DDROW,DDCOL,N = ',DDROW,DDCOL,N
              STOP
         END IF

         IF ( N.GT.DDVEC )  THEN
              WRITE (1,*) ' Dimension of vector D too small: '
              WRITE (1,*) ' mat__diagonalize_real_symmetric '
              WRITE (1,*) ' DDVEC,N = ',DDVEC,N
              WRITE (*,*) ' Dimension of vector D too small: '
              WRITE (*,*) ' mat__diagonalize_real_symmetric '
              WRITE (*,*) ' DDVEC,N = ',DDVEC,N
              STOP
         END IF
C
C
C             ...check also, if dimension of array E, which will
C                contain the offdiagonals of the tridiagonal matrix,
C                is ok.
C
C
         IF ( N.GT.LDN )  THEN
              WRITE (1,*) ' Dimension of array E too small: '
              WRITE (1,*) ' mat__diagonalize_real_symmetric '
              WRITE (1,*) ' LDN, Order of matrix N = ',LDN,N
              WRITE (*,*) ' Dimension of array E too small: '
              WRITE (*,*) ' mat__diagonalize_real_symmetric '
              WRITE (*,*) ' LDN, Order of matrix N = ',LDN,N
              STOP
         END IF
C
C
C             ...handle special case, if order of matrix is 1.
C
C
         IF ( N.EQ.1 ) THEN
              D (1) = X (1,1)
              X (1,1) = 1.D0
              RETURN
         END IF
C
C
C             ...perform Householder reduction to tridiagonal form.
C
C
         DO  150  NI = 2,N
 
             I = N+2-NI
             L = I-2
             H = 0.D0
             G = X (I,I-1)

             IF (L) 140,140,20

   20        DO  30  K = 1,L
                 H = H + X(I,K)**2
   30        CONTINUE
 
             S = H+G*G
             IF ( S.GE.TOL ) GOTO 50

   40        H = 0.D0
             GOTO 140

   50        IF (H) 140,140,60

   60        L = L+1
             F = G 
             G = DSQRT(S)

             IF (F) 75,75,70

   70        G = -G
   75        H = S-F*G
             X (I,I-1) = F-G
             F = 0.D0
 
             DO  110  J = 1,L
                 X (J,I) = X (I,J) / H
                 S = 0.D0
                 DO  80  K = 1,J
                     S = S + X(J,K)*X(I,K)
   80            CONTINUE
                 J1 = J+1
                 IF ( J1.GT.L ) GOTO 100
                 DO  90  K = J1,L
                     S = S + X(K,J)*X(I,K)
   90            CONTINUE
  100            E (J) = S/H
                 F = F + S*X(J,I)
  110        CONTINUE

             F = F/(H+H)

             DO  120  J=1,L
                 E (J) = E (J) - F*X(I,J)
  120        CONTINUE

             DO  130  J=1,L
                 F = X(I,J)
                 S = E(J)
                 DO 135 K=1,J
                    X (J,K) = X (J,K) - F*E(K) - X(I,K)*S
  135            CONTINUE
  130        CONTINUE

  140        D (I) = H
             E (I-1) = G
  150    CONTINUE
C
C
C             ...accumulation of transformation matrix and intermediate
C                D vector.
C
C
  160    D (1) = X (1,1)
         X (1,1) = 1.D0

         DO  220  I=2,N
             L = I-1

             IF ( D(I) ) 200,200,170

  170        DO  190  J=1,L
                 S = 0.D0
                 DO  180  K=1,L
                     S = S + X (I,K) * X (K,J)
  180            CONTINUE
                 DO  195  K=1,L
                     X (K,J) = X (K,J) - S * X (K,I)
  195            CONTINUE
  190        CONTINUE

  200        D (I) = X (I,I)
             X (I,I) = 1.D0
             DO  210  J=1,L
                 X (I,J) = 0.D0
                 X (J,I) = 0.D0
  210        CONTINUE  
  220     CONTINUE
C
C
C             ...QL iterations.
C
C
          B = 0.D0
          F = 0.D0
          E (N) = 0.D0

          DO  340  L=1,N
              H = EPS * ( DABS(D(L)) + DABS(E(L)) )
              IF (H.GT.B) B = H

              DO  240  J=L,N
                  IF ( DABS(E(J)).LE.B ) GOTO 250
  240         CONTINUE

  250         IF ( J.EQ.L ) THEN
                   D(L) = D(L)+F
                   GOTO 340
              END IF

  260         G = D(L)
              P = (D(L+1)-G) * 0.5D0/E(L)
              R = DSQRT (P*P+1.D0)

              IF (P) 270,280,280

  270         P = P-R
              GOTO 290
  280         P = P+R
  290         D(L) = E(L)/P
              H = G - D(L)
              K = L+1

              DO  300  I=K,N
                  D(I) = D(I) - H
  300         CONTINUE
              F = F+H
              P = D(J)
              C = 1.D0
              S = 0.D0
              J1 = J-1

              DO  330  NI=L,J1
                  I = L+J1-NI
                  G = C * E(I)
                  H = C*P
                  IF ( DABS(P).LT.DABS(E(I)) ) THEN
                       C = P/E(I)
                       R = DSQRT (C*C+1.D0)
                       E(I+1) = S*E(I)*R
                       S = 1.D0/R
                       C = C/R
                  ELSE
                       C = E(I)/P
                       R = DSQRT (C*C+1.D0)
                       E(I+1) = S*P*R
                       S = C/R
                       C = 1.D0/R
                  END IF

                  P = C*D(I) - S*G
                  D(I+1) = H + S*(C*G+S*D(I))

                  DO  335  K=1,N
                      H = X (K,I+1)
                      X (K,I+1) = X (K,I)*S + H*C
                      X (K,I) = X (K,I)*C  - H*S
  335             CONTINUE
  330         CONTINUE

              E(L) = S*P
              D(L) = C*P

              IF ( DABS(E(L)).GT.B ) GOTO 260
              D(L) = D(L)+F
  340     CONTINUE
C
C
C             ...put eigenvalues and eigenvectors in desired order.
C
C
          NI = N-1
          DO  420  I=1,NI 
              K = I
              P = D(I)
              J1 = I+1

              IF (REVERS) THEN
                  DO  360  J=J1,N 
                      IF ( D(J).LE.P ) GOTO 360 
                      K = J
                      P = D(J)
  360             CONTINUE
              ELSE
                  DO  370  J=J1,N 
                      IF ( D(J).GE.P ) GOTO 370 
                      K = J
                      P = D(J)
  370             CONTINUE
              END IF

              IF ( K.EQ.I ) GOTO 420

              D(K) = D(I)
              D(I) = P

              DO  410  J=1,N
                  P = X(J,I)
                  X(J,I) = X(J,K)
                  X(J,K) = P
  410         CONTINUE
  420     CONTINUE
C
C
C             ...finished!
C
C
          RETURN
          END
