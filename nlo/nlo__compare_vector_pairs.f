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
         SUBROUTINE  NLO__COMPARE_VECTOR_PAIRS
     +
     +                    ( DDV,DDW,DDP,
     +                      N,
     +                      ABSOLUT,
     +                      THRESH,
     +                      P,
     +                      V,W,
     +
     +                             MATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__COMPARE_VECTOR_PAIRS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : Given two vector pairs V1,V2 and W1,W2, this routine
C                examines, if the components of V1 and W1 and the
C                components of V2 and W2 match to within a threshold
C                accuracy. There is also the option to match only
C                the absolute component values. The procedure is a
C                simple exhaustive comparison between the components.
C
C
C                  Input:
C
C                    DDV,DDW   =  row dimensions of vector pairs V1,V2
C                                 and W1,W2
C                    DDP       =  dimension of index vector P
C                    N         =  # of rows of vector pairs V1,V2 and
C                                 W1,W2
C                    ABSOLUT   =  is true, if the absolute values
C                                 of the vector components are to be
C                                 taken for comparison.
C                    P         =  integer index vector that will hold
C                                 indices of vector pair V1,V2 still
C                                 to be compared with.
C                    V,W       =  N x 2 arrays corresponding to the
C                                 vector pairs V1,V2 and W1,W2.
C
C                  Output:
C
C                    MATCH     =  is true, if the components of the
C                                 vector pairs were found to match,
C                                 false otherwise.
C
C
C                Note: Both vector pairs V1,V2 and W1,W2 will NOT
C                      be destroyed!
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

         LOGICAL     ABSOLUT
         LOGICAL     EQUAL
         LOGICAL     MATCH

         INTEGER     DDV,DDW,DDP
         INTEGER     I,J,K,N
         INTEGER     LEFT

         INTEGER     P (1:DDP)

         DOUBLE PRECISION  D1,D2
         DOUBLE PRECISION  V1,V2,W1,W2
         DOUBLE PRECISION  THRESH

         DOUBLE PRECISION  V (1:DDV,1:2)
         DOUBLE PRECISION  W (1:DDW,1:2)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of vectors supplied.
C
C
         IF (N.GT.DDV .OR. N.GT.DDW) THEN
             WRITE (*,*) ' Dimensions of vector V and/or W too small! '
             WRITE (*,*) ' nlo__compare_vector_pairs '
             WRITE (*,*) ' DDV,DDW,N = ',DDV,DDW,N
             WRITE (1,*) ' Dimensions of vector V and/or W too small! '
             WRITE (1,*) ' nlo__compare_vector_pairs '
             WRITE (1,*) ' DDV,DDW,N = ',DDV,DDW,N
             STOP
         END IF

         IF (N.GT.DDP) THEN
             WRITE (*,*) ' Dimensions of index vector P too small! '
             WRITE (*,*) ' nlo__compare_vector_pairs '
             WRITE (*,*) ' DDP,N = ',DDP,N
             WRITE (1,*) ' Dimensions of index vector P too small! '
             WRITE (1,*) ' nlo__compare_vector_pairs '
             WRITE (1,*) ' DDP,N = ',DDP,N
             STOP
         END IF
C
C
C             ...proceed with exhaustive search.
C
C
         DO I = 1,N
            P (I) = I
         END DO

         LEFT = N
         MATCH = .TRUE.
C
C
C             ...based on absolute vector component values.
C
C
         IF (ABSOLUT) THEN

             DO 100 I = 1,N

                W1 = DABS (W (I,1))
                W2 = DABS (W (I,2))

                DO J = 1,LEFT
                   K = P (J)
                   V1 = DABS (V (K,1))
                   V2 = DABS (V (K,2))
                   D1 = DABS (V1 - W1)
                   D2 = DABS (V2 - W2)
                   EQUAL = D1.LT.THRESH .AND. D2.LT.THRESH
                   IF (EQUAL) THEN
                       DO K = J+1,LEFT
                          P (K-1) = P (K)
                       END DO
                       GOTO 100
                   END IF
                END DO

                MATCH = .FALSE.
                RETURN

  100        CONTINUE

         ELSE
C
C
C             ...based on true vector component values.
C
C
             DO 200 I = 1,N

                W1 = W (I,1)
                W2 = W (I,2)

                DO J = 1,LEFT
                   K = P (J)
                   V1 = V (K,1)
                   V2 = V (K,2)
                   D1 = DABS (V1 - W1)
                   D2 = DABS (V2 - W2)
                   EQUAL = D1.LT.THRESH .AND. D2.LT.THRESH
                   IF (EQUAL) THEN
                       DO K = J+1,LEFT
                          P (K-1) = P (K)
                       END DO
                       GOTO 200
                   END IF
                END DO

                MATCH = .FALSE.
                RETURN

  200        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
