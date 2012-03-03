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

         SUBROUTINE  NLO__WSW_ORTHONORMALIZE
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWS,DDCOLS,
     +                      DDVECW,
     +                      N,M,
     +                      S,W,
     +                      LOWDIN,
     +                      NOOVLP,
     +                      SAVES,
     +                      XVEC1,XVEC2,
     +                      XMAT,
     +
     +                              FAILED,
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__WSW_ORTHONORMALIZE
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine performs a weighted orthonormalization
C                on a set of M molecular orbitals, expanded by a
C                N x M coefficient matrix in a set of N atomic orbitals.
C
C                Procedure:
C
C                Given the overlap matrix S in te AO basis and the
C                M weights associated with the MO's, we first generate
C                the so called weighted overlap MO matrix
C
C                                  S(w) = WSW
C
C                where W is the diagonal weight matrix and S is the
C                MO overlap matrix. The final transformation matrix
C                O(w) for calculating the new weighted orthonormalized
C                MO's is given by
C
C                                            -1/2
C                              O(w) = W * S(w)
C
C                where the invers square root of the weighted overlap
C                MO matrix is obtained by standard methods. Note that
C                if S is positive definite, that is for all possible
C                vectors c we have
C
C                           c(T)Sc > 0       T = transpose
C
C                then also WSW must be positive definite, which can
C                be shown by considering a transformed vector c' = Wc
C                and from the fact that W is diagonal
C
C                      c'(T)Sc' = c(T)W(T)SWc = c(T)(WSW)c > 0
C
C                If LOWDIN is set true, then the routine performs only
C                a symmetric orthonormalization. In this case the
C                weight vector passed is arbitrary and all steps
C                involving the weights are skipped. Hence this would
C                correspond to a transformation matrix of the form:
C
C                                               -1/2
C                              O(w=constant) = S
C
C                If NOOVLP is set true, the routine assumes a unit
C                AO overlap matrix, with resulting computational
C                simplifications. In this case whatever contents
C                matrix S has on input will not change and the keyword
C                SAVES can be arbitrarily set.
C
C                If the routine encounters difficulties in the sense
C                of too small eigenvalues when evaluating the invers
C                square root of the (weighted) overlap matrix, the
C                routine exits with FAILED = .true.
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

         LOGICAL     FAILED
         LOGICAL     LTRG,UTRG
         LOGICAL     LOWDIN
         LOGICAL     NOOVLP
         LOGICAL     REVERS
         LOGICAL     SAVES

         INTEGER     DDROWC,DDCOLC,DDROWS,DDCOLS,DDVECW
         INTEGER     I,J,K,M,N
         INTEGER     JEND

         DOUBLE PRECISION  ZERO,ONE,SMALL
         DOUBLE PRECISION  X

         DOUBLE PRECISION  W     (1:DDVECW)
         DOUBLE PRECISION  XVEC1 (1:M)
         DOUBLE PRECISION  XVEC2 (1:M)

         DOUBLE PRECISION  C     (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  S     (1:DDROWS,1:DDCOLS)
         DOUBLE PRECISION  XMAT  (1:N,1:M)

         DATA  ONE     /1.D0/
         DATA  ZERO    /0.D0/
         DATA  SMALL   /1.0D-10/
         DATA  REVERS  /.TRUE./
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of C matrix supplied.
C
C
         IF (N.GT.DDROWC .OR. M.GT.DDCOLC) THEN
             WRITE (*,*) ' Dimensions of matrix C too small! '
             WRITE (*,*) ' nlo__wsw_orthonormalize '
             WRITE (*,*) ' DDROWC,DDCOLC,N,M = ',DDROWC,DDCOLC,N,M
             WRITE (1,*) ' Dimensions of matrix C too small! '
             WRITE (1,*) ' nlo__wsw_orthonormalize '
             WRITE (1,*) ' DDROWC,DDCOLC,N,M = ',DDROWC,DDCOLC,N,M
             STOP
         END IF
C
C
C             ...calculate lower triangle of M x M weighted overlap
C                matrix in the MO basis using the original overlap
C                matrix S in AO basis, the MO coefficient matrix C
C                in AO basis and the weight vector W. Check dimensions
C                of the overlap matrix and weight vector, if necessary.
C
C
         IF (NOOVLP) THEN
             DO 10 J = 1,M
             DO 10 I = J,M
                X = ZERO
                DO 20 K = 1,N
                   X = X + C (K,I) * C (K,J)
   20           CONTINUE
                XMAT (I,J) = X
   10        CONTINUE
         ELSE
             IF (N.GT.DDROWS .OR. N.GT.DDCOLS) THEN
                 WRITE (*,*) ' Dimensions of matrix S too small! '
                 WRITE (*,*) ' nlo__wsw_orthonormalize '
                 WRITE (*,*) ' DDROWS,DDCOLS,N = ',DDROWS,DDCOLS,N
                 WRITE (1,*) ' Dimensions of matrix S too small! '
                 WRITE (1,*) ' nlo__wsw_orthonormalize '
                 WRITE (1,*) ' DDROWS,DDCOLS,N = ',DDROWS,DDCOLS,N
                 STOP
             END IF

             LTRG  = .TRUE.
             UTRG  = .FALSE.

             CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                  ( DDROWC,DDCOLC,
     +                    DDROWS,DDCOLS,
     +                    N,M,
     +                    M,
     +                    M,N,
     +                    0,
     +                    SAVES,LTRG,UTRG,
     +                    C,S,
     +                    XVEC1,
     +
     +                             XMAT )
     +
     +
         END IF

         IF (.NOT.LOWDIN) THEN
              IF (M.GT.DDVECW) THEN
                  WRITE (*,*) ' Dimension of vector W too small! '
                  WRITE (*,*) ' nlo__wsw_orthonormalize '
                  WRITE (*,*) ' DDVECW,M = ',DDVECW,M
                  WRITE (1,*) ' Dimension of vector W too small! '
                  WRITE (1,*) ' nlo__wsw_orthonormalize '
                  WRITE (1,*) ' DDVECW,M = ',DDVECW,M
                  STOP
              END IF

              DO 100 J = 1,M
                 X = W (J)
                 DO 110 I = J,M
                    XMAT (I,J) = W (I) * X * XMAT (I,J)
  110            CONTINUE
  100         CONTINUE
         END IF
C
C
C             ...determine invers square root of weighted M x M
C                overlap matrix and calculate final transformation
C                matrix O(w).
C
C
         CALL    MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                ( N,M,M,
     +                  M,
     +                  REVERS,
     +
     +                          XVEC1,
     +                          XMAT )
     +
     +
         FAILED = .FALSE.

         DO 200 I = 1,M
            X = XVEC1 (I)
            IF ( X.LT.SMALL ) THEN
                 WRITE (*,*) ' M x M S(w) matrix not +ve definite! '
                 WRITE (*,*) ' EIGVAL,SMALL = ',X,SMALL
                 WRITE (*,*) ' nlo__wsw_orthonormalize '
                 WRITE (1,*) ' M x M S(w) matrix not +ve definite! '
                 WRITE (1,*) ' EIGVAL,SMALL = ',X,SMALL
                 WRITE (1,*) ' nlo__wsw_orthonormalize '
                 FAILED = .TRUE.
                 RETURN
            END IF
            XVEC1 (I) = ONE / DSQRT (X)
  200    CONTINUE

         CALL    MAT__C_EQ_C_TRANSPOSED_FLOAT
     +
     +                ( N,M,
     +                  M,M,
     +                  XMAT )
     +
     +
         DO 210 J = 1,M
            DO 220 I = J,M
               X = ZERO
               DO 230 K = 1,M
                  X = X + XVEC1 (K) * XMAT (K,I) * XMAT (K,J)
  230          CONTINUE
               XVEC2 (I) = X
  220       CONTINUE
            DO 240 I = J,M
               XMAT (I,J) = XVEC2 (I)
  240       CONTINUE
  210    CONTINUE

         IF (LOWDIN) THEN
             DO 250 I = 1,M
                JEND = I - 1
                DO 260 J = 1,JEND
                   XMAT (J,I) = XMAT (I,J)
  260           CONTINUE
  250        CONTINUE
         ELSE
             DO 270 I = 1,M
                X = W (I)
                JEND = I - 1
                DO 280 J = 1,JEND
                   XMAT (J,I) = W (J) * XMAT (I,J)
                   XMAT (I,J) =     X * XMAT (I,J)
  280           CONTINUE
                XMAT (I,I) = X * XMAT (I,I)
  270        CONTINUE
         END IF
C
C
C             ...update N x M pre-NAO coefficient matrix with
C                M x M transformation matrix O(w) to obtain the
C                M NAO's.
C
C
         CALL    MAT__C_EQ_C_TIMES_A_FLOAT
     +
     +                ( DDROWC,DDCOLC,
     +                  N,M,
     +                  M,M,
     +                  N,M,M,
     +                  XVEC1,XVEC2,
     +                  XMAT,
     +
     +                          C )
     +
C
C
C             ...ready!
C
C
         RETURN
         END
