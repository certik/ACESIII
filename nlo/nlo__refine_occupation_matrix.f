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
         SUBROUTINE  NLO__REFINE_OCCUPATION_MATRIX
     +
     +                    ( DDROWP,DDCOLP,
     +                      N,
     +                      MAXOCC,
     +                      QSYMACC,
     +
     +                              P )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__REFINE_OCCUPATION_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine tries to remove certain elements from
C                the occupation matrix which might have crept in due
C                to rounding errors. These rounding errors show up
C                as offdiagonal noise and can be eliminated for certain
C                rows and columns.
C
C                The N x N occupation matrix P obeys the idempotency
C                condition:
C
C                                 P = P * P
C
C                For any diagonal element in P we have thus:
C
C                                                N
C                    P (I,I) - P (I,I)**2  =    sum   P (I,K)**2
C                                             i neq k
C
C                Consider first a small diagonal element P (I,I) = E,
C                where E is so small that E*E can be neglected. Then
C                we have:
C                                         N
C                                 E =    sum   P (I,K)**2
C                                      i neq k
C
C                Now let us say that the maximum absolute value of all
C                the P (I,K)'s is less than T / sqrt (N-1):
C
C                       |P (I,K)|     <  T / sqrt (N-1)
C                                max
C
C                where T will be determined shortly. Then we have
C                the following inequality:
C
C                            N
C                           sum   P (I,K)**2   <  T*T
C                         i neq k
C
C                Thus if we set T = sqrt (E), the diagonal element
C                cannot exceed the value of E. Setting thus T = 1.D-10
C                and if the maximum |P (I,K)| value is below the value
C                of T / sqrt (N-1), the diagonal element can only be
C                < 1.D-20, which is below the accuracy obtained by
C                the matrix diagonalization procedures. If this
C                condition on the offdiagonal elements is found we
C                can safely set all the offdiagonals equal to zero.
C                The diagonal element will be set equal to MAXOCC or
C                zero. The diagonal MAXOCC case can be analyzed the
C                same way as the small diagonal E case, since, if
C                P (I,I) = MAXOCC - E, then we have:
C
C                       P (I,I) - P (I,I)**2  =  E - E*E
C
C                and the same reasoning as above follows through.
C
C
C                  Input:
C
C                    DDROWP,DDCOLP  =  declared dimensions of matrix P
C                    N              =  current dimension of occupation
C                                      matrix P.
C                    MAXOCC         =  maximum orbital occupancy number
C                                      (can be only 1 or 2).
C                    QSYMACC      =  symmetry accuracy for orbital
C                                    interaction order values (sum of
C                                    squares of offdiagonal occupation
C                                    matrix elements for one orbital)
C                    P              =  original occupation matrix.
C                                      Only lower triangle needs to be 
C                                      supplied.
C
C                  Output:
C
C                    P              =  lower triangle of refined
C                                      occupation matrix.
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

C         LOGICAL     WHIGH,WZERO

         INTEGER     DDROWP,DDCOLP
         INTEGER     I,J,N

         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  Q,QSYMACC
C         DOUBLE PRECISION  PMAX
C         DOUBLE PRECISION  W,WDIFF
         DOUBLE PRECISION  T
C         DOUBLE PRECISION  X
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  P (1:DDROWP,1:DDCOLP)

         DATA  T     /1.D-10/
         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of P matrix supplied.
C
C
         IF (N.GT.DDROWP .OR. N.GT.DDCOLP) THEN
             WRITE (*,*) ' Dimensions of matrix P too small! '
             WRITE (*,*) ' nlo__refine_occupation_matrix '
             WRITE (*,*) ' DDROWP,DDCOLP,N = ',DDROWP,DDCOLP,N
             WRITE (1,*) ' Dimensions of matrix P too small! '
             WRITE (1,*) ' nlo__refine_occupation_matrix '
             WRITE (1,*) ' DDROWP,DDCOLP,N = ',DDROWP,DDCOLP,N
             STOP
         END IF
C
C
C             ...form the interaction orders Q and zero out the
C                columns and rows of P which correspond to Q values
C                below the symmetry accuracy.
C
C
         DO J = 1,N

            Q = ZERO
            DO I = J+1,N
               Q = Q + P (I,J) ** 2
            END DO
            DO I = J-1,1,-1
               Q = Q + P (J,I) ** 2
            END DO
            Q = DSQRT (Q)

            IF (Q.LT.QSYMACC) THEN
                DO I = J+1,N
                   P (I,J) = ZERO
                END DO
                DO I = J-1,1,-1
                   P (J,I) = ZERO
                END DO
            END IF

         END DO

C         GOTO 1234
C
C
C             ...loop over all diagonals (weights) of P.
C
C
C         X = T / DSQRT (DFLOAT (N-1))
C
C         DO J = 1,N
C
C            PMAX = ZERO
C            DO I = J+1,N
C               PMAX = MAX (PMAX,DABS (P (I,J)))
C            END DO
C            DO I = J-1,1,-1
C               PMAX = MAX (PMAX,DABS (P (J,I)))
C            END DO
C
C            IF (PMAX.LT.X) THEN
C
C                DO I = J+1,N
C                   P (I,J) = ZERO
C                END DO
C                DO I = J-1,1,-1
C                   P (J,I) = ZERO
C                END DO
C
C                W = DABS (P (J,J))
C
C                WZERO = W .LT. T
C                WDIFF = MAXOCC - W
C                WHIGH = DABS (WDIFF) .LT. T
C
C                IF (WZERO.AND.WHIGH) THEN
C                    WRITE (*,*) ' Problems refining P matrix! '
C                    WRITE (*,*) ' WZERO,WHIGH = ',WZERO,WHIGH
C                    WRITE (*,*) ' nlo__refine_occupation_matrix '
C                    WRITE (1,*) ' Problems refining P matrix! '
C                    WRITE (1,*) ' WZERO,WHIGH = ',WZERO,WHIGH
C                    WRITE (1,*) ' nlo__refine_occupation_matrix '
C                    STOP
C                ELSE IF (WZERO) THEN
C                    P (J,J) = ZERO
C                ELSE IF (WHIGH) THEN
C                    P (J,J) = MAXOCC
C                ELSE
C                    WRITE (*,*) ' Problems refining P matrix! '
C                    WRITE (*,*) ' WZERO,WHIGH = ',WZERO,WHIGH
C                    WRITE (*,*) ' nlo__refine_occupation_matrix '
C                    WRITE (1,*) ' Problems refining P matrix! '
C                    WRITE (1,*) ' WZERO,WHIGH = ',WZERO,WHIGH
C                    WRITE (1,*) ' nlo__refine_occupation_matrix '
C                    STOP
C                END IF
C
C            END IF
C
C         END DO
C
C
C             ...ready!
C
C
C 1234    RETURN
C
         RETURN
         END
