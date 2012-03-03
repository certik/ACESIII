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
         SUBROUTINE  MAT__BLOCKDIAG_REAL_SYMMETRIC
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWY,DDCOLY,
     +                      DDVECX,DDVECY,
     +                      DDVECD,
     +                      N,
     +                      OFFZERO,
     +                      REVERS,
     +                      IX,IY,
     +                      X,
     +
     +                              D,
     +                              Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__BLOCKDIAG_REAL_SYMMETRIC
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This routine performs a block-diagonalization on a real
C                symmetric matrix X. The criterion for establishing
C                the different blocks is OFFZERO, which is the smallest
C                absolute value below which an offdiagonal matrix
C                element of X is considered to be zero. The individual
C                subblocks of X are set up individually and diagonalized
C                separately. The final eigenvectors of X are assembled
C                with zeros in the appropriate places corresponding
C                to noninteracting matrix elements of X.
C
C
C                  Input:
C
C                   DDROWX = row dimension for matrix X to be block-
C                            diagonalized
C                   DDCOLX = column dimension for matrix X to be block-
C                            diagonalized
C                   DDROWY = row dimension for matrix Y to hold the
C                            individual subblocks
C                   DDCOLY = column dimension for matrix Y to hold the
C                            individual subblocks
C                   DDVECX = dimension for index vector IX, which will
C                            contain the info to which subblock # each
C                            row and column of X will belong to
C                   DDVECY = dimension for index vector IY, which will
C                            hold the row and column indices of X for
C                            each subblock
C                   DDVECD = dimension for vector D containing the
C                            eigenvalues
C                        N = order of matrix X to be block-diagonalized
C                  OFFZERO = offdiagonal zero matrix element threshold
C                   REVERS = if false, eigenvalues and associated
C                            eigenvectors are in ascending sequence,
C                            otherwise this ordering is reversed.
C                    IX,IY = index vectors
C                        X = matrix to be block-diagonalized (only
C                            lower triangle needs to be supplied)
C
C                  Output:
C
C                        D = computed eigenvalues of X, in specified
C                            order
C                        Y = matrix of corresponding eigenvectors,
C                            columnwise.
C
C                Note, that the original X is NOT! destroyed during the
C                process!!!
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

         LOGICAL     ADD
         LOGICAL     CASE1,CASE2
         LOGICAL     CHANGE
         LOGICAL     REVERS

         INTEGER     COL,ROW
         INTEGER     DDROWX,DDCOLX,DDROWY,DDCOLY
         INTEGER     DDVECX,DDVECY,DDVECD
         INTEGER     I,J,K,L,M,N
         INTEGER     LAST
         INTEGER     MSUB,NSUB
         INTEGER     NSWEEP
         INTEGER     SWEEP

         INTEGER     IX (1:DDVECX)
         INTEGER     IY (1:DDVECY)

         DOUBLE PRECISION  E
         DOUBLE PRECISION  ZERO,OFFZERO

         DOUBLE PRECISION  D (1:DDVECD)

         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)
         DOUBLE PRECISION  Y (1:DDROWY,1:DDCOLY)

         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...check passed dimension of matrices and vectors.
C
C
         IF (N.GT.DDROWX .OR. N.GT.DDCOLX)  THEN
             WRITE (*,*) ' Dimension of matrix X too small: '
             WRITE (*,*) ' mat__blockdiag_real_symmetric '
             WRITE (*,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             WRITE (1,*) ' Dimension of matrix X too small: '
             WRITE (1,*) ' mat__blockdiag_real_symmetric '
             WRITE (1,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             STOP
         END IF

         IF (N.GT.DDROWY .OR. N.GT.DDCOLY)  THEN
             WRITE (*,*) ' Dimension of matrix Y too small: '
             WRITE (*,*) ' mat__blockdiag_real_symmetric '
             WRITE (*,*) ' DDROWY,DDCOLY,N = ',DDROWY,DDCOLY,N
             WRITE (1,*) ' Dimension of matrix Y too small: '
             WRITE (1,*) ' mat__blockdiag_real_symmetric '
             WRITE (1,*) ' DDROWY,DDCOLY,N = ',DDROWY,DDCOLY,N
             STOP
         END IF

         IF (N.GT.DDVECX)  THEN
             WRITE (*,*) ' Dimension of vector IX too small: '
             WRITE (*,*) ' mat__blockdiag_real_symmetric '
             WRITE (*,*) ' DDVECX,N = ',DDVECX,N
             WRITE (1,*) ' Dimension of vector IX too small: '
             WRITE (1,*) ' mat__blockdiag_real_symmetric '
             WRITE (1,*) ' DDVECX,N = ',DDVECX,N
             STOP
         END IF

         IF (N.GT.DDVECY)  THEN
             WRITE (*,*) ' Dimension of vector IY too small: '
             WRITE (*,*) ' mat__blockdiag_real_symmetric '
             WRITE (*,*) ' DDVECY,N = ',DDVECY,N
             WRITE (1,*) ' Dimension of vector IY too small: '
             WRITE (1,*) ' mat__blockdiag_real_symmetric '
             WRITE (1,*) ' DDVECY,N = ',DDVECY,N
             STOP
         END IF

         IF (N.GT.DDVECD)  THEN
             WRITE (*,*) ' Dimension of vector D too small: '
             WRITE (*,*) ' mat__blockdiag_real_symmetric '
             WRITE (*,*) ' DDVECD,N = ',DDVECD,N
             WRITE (1,*) ' Dimension of vector D too small: '
             WRITE (1,*) ' mat__blockdiag_real_symmetric '
             WRITE (1,*) ' DDVECD,N = ',DDVECD,N
             STOP
         END IF
C
C
C             ...detect the presence of submatrices. MSUB will contain
C                the total # of such submatrices and IX will have the
C                submatrix number on those rows of matrix X belonging
C                together.
C
C
         DO ROW = 1,N
            IX (ROW) = 0
         END DO

         MSUB = 0

         DO ROW = 1,N
            IF (IX (ROW).EQ.0) THEN

                MSUB = MSUB + 1
                IX (ROW) = MSUB

                COL = ROW
                DO I = ROW+1,N
                   CASE1 =          IX (I) .EQ. 0
                   CASE2 = ABS (X (I,COL)) .GE. OFFZERO
                   IF (CASE1.AND.CASE2) THEN
                       IX (I) = MSUB
                   END IF
                END DO

                CHANGE = .TRUE.
                NSWEEP = N - ROW

                DO SWEEP = 1,NSWEEP
                   IF (CHANGE) THEN
                       CHANGE = .FALSE.
                       DO COL = ROW+1,N
                          ADD = .FALSE.
                          DO I = COL,N
                             CASE1 =          IX (I) .EQ. MSUB
                             CASE2 = ABS (X (I,COL)) .GE. OFFZERO
                             ADD = ADD .OR. (CASE1.AND.CASE2)
                          END DO
                          IF (ADD) THEN
                              DO I = COL,N
                                 IF (ABS (X (I,COL)).GE.OFFZERO) THEN
                                     IX (I) = MSUB
                                     CHANGE = .TRUE.
                                 END IF
                              END DO
                          END IF
                       END DO
                   END IF
                END DO

            END IF
         END DO

         DO ROW = 1,N
            IF (IX (ROW).EQ.0) THEN
                WRITE (*,*) ' Problems in finding submatrices! '
                WRITE (*,*) ' mat__blockdiag_real_symmetric '
                WRITE (*,*) ' ROW,IX (ROW) = ',ROW,IX (ROW)
                WRITE (1,*) ' Problems in finding submatrices! '
                WRITE (1,*) ' mat__blockdiag_real_symmetric '
                WRITE (1,*) ' ROW,IX (ROW) = ',ROW,IX (ROW)
                STOP
            END IF
         END DO
C
C
C             ...if only one submatrix, perform straightforward overall
C                diagonalization.
C
C
         IF (MSUB.EQ.1) THEN

             CALL  MAT__C_EQ_A_FLOAT
     +
     +                  ( DDROWX,DDCOLX,
     +                    DDROWY,DDCOLY,
     +                    N,N,
     +                    X,
     +
     +                            Y )
     +
     +
             CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                  ( DDROWY,DDCOLY,DDVECD,
     +                    N,
     +                    REVERS,
     +
     +                            D,
     +                            Y )
     +
     +
             RETURN
         END IF
C
C
C             ...loop over all submatrices, assemble them, diagonalize
C                them and form the overall eigenvectors of X.
C
C
         LAST = 0

         DO M = 1,MSUB

            NSUB = 0
            DO I = 1,N
               IF (IX (I).EQ.M) THEN
                   NSUB = NSUB + 1
                   IY (NSUB) = I
               END IF
            END DO

            DO J = 1,NSUB
               COL = IY (J)
               DO I = J,NSUB
                  ROW = IY (I)
                  Y (I,LAST+J) = X (ROW,COL)
               END DO
            END DO

            CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                 ( DDROWY,NSUB,NSUB,
     +                   NSUB,
     +                   REVERS,
     +
     +                           D (LAST+1),
     +                           Y (1,LAST+1) )
     +
     +
            DO J = 1,NSUB
               L = N
               COL = LAST + J
               DO I = NSUB,1,-1
                  ROW = IY (I)
                  DO K = L,ROW+1,-1
                     Y (K,COL) = ZERO
                  END DO
                  Y (ROW,COL) = Y (I,COL)
                  L = ROW - 1
               END DO
               DO K = L,1,-1
                  Y (K,COL) = ZERO
               END DO
            END DO

            LAST = LAST + NSUB

         END DO
C
C
C             ...put eigenvalues and eigenvectors of matrix X in
C                desired order.
C
C
         DO I = 1,N-1
            K = I
            E = D (I)

            IF (REVERS) THEN
                DO J = I+1,N
                   IF (D (J).GT.E) THEN
                       K = J
                       E = D (J)
                   END IF
                END DO
            ELSE
                DO J = I+1,N
                   IF (D (J).LT.E) THEN
                       K = J
                       E = D (J)
                   END IF
                END DO
            END IF

            IF (K.NE.I) THEN
                D (K) = D (I)
                D (I) = E
                DO J = 1,N
                   E = Y (J,I)
                   Y (J,I) = Y (J,K)
                   Y (J,K) = E
                END DO
            END IF

         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
