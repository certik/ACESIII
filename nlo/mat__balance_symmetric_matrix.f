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
         SUBROUTINE  MAT__BALANCE_SYMMETRIC_MATRIX
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDVECW,
     +                      N,
     +                      THRESH,
     +                      W,
     +
     +                              X )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__BALANCE_SYMMETRIC_MATRIX
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : Given a symmetric matrix X, this routine balances its
C                elements to within a threshold THRESH. By balancing
C                we mean setting all elements strictly equal that are
C                equal to within the threshold value.
C
C
C                Procedure:
C                ----------
C
C                  Loop over all matrix elements (lower triangle)
C                  column by column
C
C                      If element X(i,j) was not balanced:
C
C                    a)   Set SUM = X(i,j)
C                               M = 1
C
C                         Loop over all the rest of the X matrix
C
C                             If ABS (X(k,l)-X(i,j)) < THRESH, then
C                                SUM = SUM + X(k,l)
C                                  M = M + 1
C                             end if
C
C                         end do
C
C                         Form average XAV = SUM / M
C
C                         Set X(i,j) = XAV and go back to step a)
C
C                         Keep doing until value of M is stable, that
C                         is equal between two runs.
C
C                         If M = stable, rerun from step a) but this
C                         time place XAV into X(i,j) and all
C                         relevant X(k,l) and mark them as balanced.
C
C                      end if
C
C                  end do
C
C
C                  Input:
C
C                   DDROWX = row dimension for matrix X
C                   DDCOLX = column dimension for matrix X
C                   DDROWW = row dimension for working matrix W
C                   DDCOLW = column dimension for working matrix W
C                   DDVECD = dimension for vector D
C                        N = order of the matrix X to be balanced
C                        D = vector D that will contain the eigenvalues
C                            of X and later the updated transfomed
C                            diagonals
C                        X = matrix to be diagonally equalized
C
C                  Output:
C
C                        X = transformation matrix
C
C                Only the lower triangle of X needs to be supplied.
C
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   BALANCE
         LOGICAL   STABLE

         INTEGER   DDROWX,DDCOLX,DDVECW
         INTEGER   I,J,K,L,M,N
         INTEGER   ITER
         INTEGER   MOLD

         DOUBLE PRECISION   SUM
         DOUBLE PRECISION   THRESH
         DOUBLE PRECISION   XAV
         DOUBLE PRECISION   ZERO,ONE

         DOUBLE PRECISION   W (1:DDVECW)

         DOUBLE PRECISION   X (1:DDROWX,1:DDCOLX)

         PARAMETER   (ZERO = 0.D0)
         PARAMETER   (ONE  = 1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check passed dimensions.
C
C
         IF (N.GT.DDROWX .OR. N.GT.DDCOLX) THEN
             WRITE (*,*) ' Dimension of matrix X too small: '
             WRITE (*,*) ' mat__balance_symmetric_matrix '
             WRITE (*,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             WRITE (1,*) ' Dimension of matrix X too small: '
             WRITE (1,*) ' mat__balance_symmetric_matrix '
             WRITE (1,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             STOP
         END IF

         IF (N.GT.DDVECW)  THEN
             WRITE (1,*) ' Dimension of vector W too small: '
             WRITE (1,*) ' mat__balance_symmetric_matrix '
             WRITE (1,*) ' DDVECW,N = ',DDVECW,N
             WRITE (*,*) ' Dimension of vector W too small: '
             WRITE (*,*) ' mat__balance_symmetric_matrix '
             WRITE (*,*) ' DDVECW,N = ',DDVECW,N
             STOP
         END IF
C
C
C             ...set upper triangle of X and W vector to zero,
C                to mark all elements of lower triangle of X as
C                unbalanced.
C
C
         DO J = 1,N
            W (J) = ZERO
            DO I = 1,J-1
               X (I,J) = ZERO
            END DO
         END DO
C
C
C             ...outer loop over all X matrix elements columnwise.
C
C
         DO J = 1,N
         DO I = J,N

            IF (I.NE.J) THEN
                BALANCE = X (J,I) .EQ. ZERO
            ELSE
                BALANCE = W (I) .EQ. ZERO
            END IF

            IF (BALANCE) THEN
                XAV = X (I,J)
                MOLD = 1
                ITER = 0

 1000           IF (ITER.GT.10) THEN
                    WRITE (*,*) ' # of balancing iterations > 10! '
                    WRITE (*,*) ' mat__balance_symmetric_matrix '
                    WRITE (1,*) ' # of balancing iterations > 10! '
                    WRITE (1,*) ' mat__balance_symmetric_matrix '
                    STOP
                END IF

                M = 1
                SUM = XAV
                DO K = I+1,N
                   IF (ABS (X (K,J) - XAV).LT.THRESH) THEN
                       M = M + 1
                       SUM = SUM + X (K,J)
                   END IF
                END DO
                DO L = J+1,N
                DO K = L,N
                   IF (ABS (X (K,J) - XAV).LT.THRESH) THEN
                       M = M + 1
                       SUM = SUM + X (K,J)
                   END IF
                END DO
                END DO

                XAV = SUM / DFLOAT (M)
                STABLE = M.EQ.MOLD

                IF (STABLE) THEN
                    X (I,J) = XAV
                    IF (I.NE.J) THEN
                        X (J,I) = ONE
                    ELSE
                        W (I) = ONE
                    END IF
                    DO K = I+1,N
                       IF (ABS (X (K,J) - XAV).LT.THRESH) THEN
                           X (K,J) = XAV
                           X (J,K) = ONE
                       END IF
                    END DO
                    DO L = J+1,N
                       IF (ABS (X (L,L) - XAV).LT.THRESH) THEN
                           X (L,L) = XAV
                           W (L) = ONE
                       END IF
                       DO K = L+1,N
                          IF (ABS (X (K,L) - XAV).LT.THRESH) THEN
                              X (K,L) = XAV
                              X (L,K) = ONE
                          END IF
                       END DO
                    END DO
                ELSE
                    MOLD = M
                    ITER = ITER + 1
                    GOTO 1000
                END IF

            END IF

         END DO
         END DO
C
C
C             ...finished!
C
C
         RETURN
         END
