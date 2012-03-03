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
         SUBROUTINE  NLO__FIND_AXIAL_ROTATION_MATRIX
     +
     +                    ( DDROW,DDCOL,
     +                      NHYB,
     +                      P,
     +
     +                              ROT )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FIND_AXIAL_ROTATION_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine calculates the rotation matrix of
C                NHYB identical hybrids which are arranged such that
C                they are equally spaced in a plane in space.
C                The rotation is against a background interaction
C                of NHYB interaction points P on one of the hybrids.
C                The goal is to maximize the interaction of the first
C                interaction point P(1).
C
C                Example:
C                
C                Consider the three sp2 hybrids. Let the absolute values
C                of the 3 interaction points on one of the sp2 hybrids
C                (denoted by o's below) be P1,P2,P3:
C
C
C                                      o
C                                      o    .P3
C                                      o   .
C                                      o  .
C                                      o .
C                                      o.    <---- rotate (o) by theta
C                                     * *
C                           P2      *     *
C                                 *         *
C                               *             *
C                             *           P1    *
C
C
C
C                          to get (P1,P2,P3 stay fixed)
C
C
C
C                                      *
C                                      P3
C                                      *
C                                      *
C                                      *
C                                      *
C                                     * o
C                                   *     o
C                                 *         o
C                               P2            P1
C                             *                 o
C
C
C
C                Let us first find the general form of such rotation
C                matrix. The three sp2 hybrids are formed by one s and
C                the two px and py functions. A general rotation T of
C                these functions (in basis order s,px,py) in the xy
C                plane is given by (c,s = cos,sin theta):
C
C
C                                     1   0   0
C
C                              T  =   0   c   s
C
C                                     0  -s   c
C
C                since the s function does not change upon rotation.
C                Also we have the defining matrix U of the three sp2
C                hybrids:
C
C
C                                      e   e   e
C
C                              U  =   2f  -f  -f
C
C                                      0   g  -g
C
C                where e = 1/sqrt(3), f = 1/sqrt(6) and c = 1/sqrt(2).
C                The general rotation matrix R of the sp2 hybrids is
C                then:
C
C                                     -1
C                                R = U   * T * U
C
C
C                which after inserting the above values becomes:
C
C
C                     1  1  1           2 -1 -1                 0  1 -1
C               1/3   1  1  1  +  c/3  -1  2 -1  +  s/sqrt(3)  -1  0  1
C                     1  1  1          -1 -1  2                 1 -1  0
C
C
C                Forming the product R*P for the new rotated P1'
C                element (the one we intend to maximize):
C
C
C                                               P1
C
C                                               P2
C
C                                               P3
C
C                              R11  R12  R13    P1'
C
C                              R11  R12  R13    ...
C
C                              R11  R12  R13    ...
C
C                we get
C
C
C            P1' = 1/3 (P1+P2+P3) + c/3 (2P1-P2-P3) + s/sqrt(3) (P2-P3)
C
C
C                which, when maximizing with respect to c under the
C                constrain c*c+s*s=1, gives for c:
C
C
C
C                                 /      (2P1 - P2 - P3)**2
C                 c(max) = +- \  / -----------------------------------
C                              \/ 3(P2 - P3)**2  +  (2P1 - P2 - P3)**2
C
C
C                 s(max) = +- sqrt (1 - c(max)*c(max)
C
C
C                The signs to be chosen are the following:
C
C                   For c(max) : +ve  if 2P1 >= P2 + P3
C                                -ve  otherwise
C
C                   For s(max) : +ve  if P2 >= P3
C                                -ve  otherwise
C
C
C                  Input:
C
C                    DDROW        =  declared row dimension of rotation
C                                    matrix
C                    DDCOL        =  declared column dimension of
C                                    rotation matrix
C                    NHYB         =  # of degenerate hybrids.
C                    P            =  the NHYB interaction points.
C
C
C                  Output:
C
C                    ROT          =  the NHYB x NHYB rotation matrix.
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

         INTEGER     DDROW,DDCOL
         INTEGER     NHYB

         DOUBLE PRECISION  C,S
         DOUBLE PRECISION  COSINE,SINE
         DOUBLE PRECISION  P1,P2,P3
         DOUBLE PRECISION  P1MAX,P1VAL
         DOUBLE PRECISION  X,Y,Z
         DOUBLE PRECISION  ZERO,THIRD,SQR3TH,ONE

         DOUBLE PRECISION  P (1:NHYB)

         DOUBLE PRECISION  ROT (1:DDROW,1:DDCOL)

         PARAMETER  (ZERO    = 0.D0                                   )
         PARAMETER  (THIRD   = 1.D0 / 3.D0                            )
         PARAMETER  (SQR3TH  = 0.57735026918962576450914878050195750D0)
         PARAMETER  (ONE     = 1.D0                                   )
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed.
C
C
         IF (NHYB.EQ.3) THEN

             P1MAX = ZERO

             P1 = P (1)
             P2 = P (2)
             P3 = P (3)

             X = P1 + P1 - P2 - P3
             Y = P2 - P3
             X = X * X
             Y = Y * Y

             C = DSQRT (X / (Y+Y+Y+X))
             S = DSQRT (ONE - C*C)
             
             X = THIRD * (ONE + C + C)
             Y = THIRD * (ONE - C) + SQR3TH * S
             Z = THIRD * (ONE - C) - SQR3TH * S

             P1VAL = ABS (P1 * X + P2 * Z + P3 * Y)

             IF (P1VAL.GT.P1MAX) THEN
                 COSINE = C
                 SINE = S
                 P1MAX = P1VAL
             END IF

             C = - C

             X = THIRD * (ONE + C + C)
             Y = THIRD * (ONE - C) + SQR3TH * S
             Z = THIRD * (ONE - C) - SQR3TH * S

             P1VAL = ABS (P1 * X + P2 * Z + P3 * Y)

             IF (P1VAL.GT.P1MAX) THEN
                 COSINE = C
                 SINE = S
                 P1MAX = P1VAL
             END IF

             C = - C
             S = - S

             X = THIRD * (ONE + C + C)
             Y = THIRD * (ONE - C) + SQR3TH * S
             Z = THIRD * (ONE - C) - SQR3TH * S

             P1VAL = ABS (P1 * X + P2 * Z + P3 * Y)

             IF (P1VAL.GT.P1MAX) THEN
                 COSINE = C
                 SINE = S
                 P1MAX = P1VAL
             END IF

             C = - C

             X = THIRD * (ONE + C + C)
             Y = THIRD * (ONE - C) + SQR3TH * S
             Z = THIRD * (ONE - C) - SQR3TH * S

             P1VAL = ABS (P1 * X + P2 * Z + P3 * Y)

             IF (P1VAL.GT.P1MAX) THEN
                 COSINE = C
                 SINE = S
                 P1MAX = P1VAL
             END IF

             X = THIRD * (ONE + COSINE + COSINE)
             Y = THIRD * (ONE - COSINE) + SQR3TH * SINE
             Z = THIRD * (ONE - COSINE) - SQR3TH * SINE

             ROT (1,1) = X
             ROT (2,1) = Z
             ROT (3,1) = Y
             ROT (1,2) = Y
             ROT (2,2) = X
             ROT (3,2) = Z
             ROT (1,3) = Z
             ROT (2,3) = Y
             ROT (3,3) = X

         ELSE

             WRITE (*,*) ' Cannot rotate > 3 axial hybrids! '
             WRITE (*,*) ' NHYB = ',NHYB
             WRITE (*,*) ' nlo__find_axial_rotation_matrix '
             WRITE (1,*) ' Cannot rotate > 3 axial hybrids! '
             WRITE (1,*) ' NHYB = ',NHYB
             WRITE (1,*) ' nlo__find_axial_rotation_matrix '
             STOP

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
