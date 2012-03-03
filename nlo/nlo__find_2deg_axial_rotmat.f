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
         SUBROUTINE  NLO__FIND_2DEG_AXIAL_ROTMAT
     +
     +                    ( P1,P2,
     +                      MAXIMZE,
     +
     +                              R11,R12,
     +                              R21,R22 )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FIND_2DEG_AXIAL_ROTMAT
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine calculates the rotation matrix of a
C                pair of degenerate NBOs such that the interaction P1
C                between the first NBO of the pair and a background
C                interaction point is maximized (largest positive value)
C                or minimized (zero value).
C
C                Procedure:
C                
C                Consider the two functions px,py. Then the interaction
C                part of the interaction matrix corresponding to the
C                first interaction point is P1 and P2. A rotation matrix
C                will bring the first interaction point to a new value
C                P1':
C
C
C                                               cos(t)   sin(t)
C
C                                              -sin(t)   cos(t)
C
C
C                                   P1    P2     P1'      ...
C
C                                  ...   ...    ...       ...
C
C                                  ...   ...    ...       ...
C
C
C                Maximizing P1' with respect to cos(t) gives the four
C                possible solutions:
C
C
C                            cos(t) = +- sqrt (P1*P1/(P1*P1 + P2*P2))
C
C                            sin(t) = +- sqrt (1 - cos(t)*cos(t))
C
C
C                of which the one leading to the maximum positive value
C                of P1' is chosen. On the other hand, zeroing P1'
C                leads to the four possible solutions:
C
C
C                            cos(t) = +- sqrt (P2*P2/(P1*P1 + P2*P2))
C
C                            sin(t) = +- sqrt (1 - cos(t)*cos(t))
C
C
C                Note, that the minimizing solution is just the cos->sin
C                and sin->cos switched maximizing solution.
C
C
C                  Input:
C
C                    P1,P2        =  the two interaction points
C                                    corresponding to the pair of
C                                    degenerate NBOs.
C                    MAXIMZE      =  is true, if the maximizing
C                                    solution is wanted, and false,
C                                    if the minimizing one is needed.
C
C
C                  Output:
C
C                    Rxy          =  the 2 x 2 rotation matrix
C                                    components (x,y=1,2).
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

         LOGICAL  MAXIMZE

         DOUBLE PRECISION  C,S
         DOUBLE PRECISION  COSINE,SINE
         DOUBLE PRECISION  P1,P2,P1NEW
         DOUBLE PRECISION  R11,R12,R21,R22
         DOUBLE PRECISION  X
         DOUBLE PRECISION  ZERO,ONE

         PARAMETER  (ZERO = 0.D0)
         PARAMETER  (ONE  = 1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed.
C
C
         IF (MAXIMZE) THEN
             X = P2 / P1
             C = DSQRT (ONE / (ONE + X*X))
             S = DSQRT (ONE - C*C)
             COSINE = DSIGN (C,P1)
             SINE = DSIGN (S,-P2)
         ELSE
             X = P1 / P2
             C = DSQRT (ONE / (ONE + X*X))
             S = DSQRT (ONE - C*C)
             COSINE = DSIGN (C,P1)
             SINE = DSIGN (S,P2)
             P1NEW = COSINE * P1 - SINE * P2
             WRITE (*,*) ' P1NEW = ',P1NEW
         END IF

         R11 = COSINE
         R12 = SINE
         R21 = - SINE
         R22 = COSINE
C
C
C             ...ready!
C
C
         RETURN
         END
