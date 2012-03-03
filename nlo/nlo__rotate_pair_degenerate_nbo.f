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
         SUBROUTINE  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                    ( NBAS,
     +                      P,
     +                      XVEC,
     +                      XMAT,
     +
     +                            C1,S1,
     +                            C2,S2,
     +                            X1,X2 )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ROTATE_PAIR_DEGENERATE_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine rotates a pair of two-fold weight
C                degenerate atomic NBOs, such that their lobes match
C                and form rectangular angles between them:
C
C
C
C                                                  x
C                       x  *                       *
C                        x *    x                  x
C                         x* x                     *
C                    ******x******   ---->   x**x**x**x**x
C                        x *x                      *
C                     x    * x                     x
C                          *  x                    *
C                                                  x
C
C
C
C                We look for a rotation matrix R of the form:
C
C
C
C                                    -----------------
C                                   | C1  S1 |        |
C                                   |        |    0   |
C                                   |-S1  C1 |        |
C                          R  =      -----------------
C                                   |        | C2  S2 |
C                                   |   0    |        |
C                                   |        |-S2  C2 |
C                                    -----------------
C
C
C
C                which brings the original occupation matrix P
C                corresponding to the pair (w is the weight of the
C                two-fold weight degenerate atomic NBOs):
C
C
C
C                                        ------
C                                       | w  0 |
C                                       |      |
C                                       | 0  w |
C                              P  =      -------------
C                                       | p1 p2| w  0 |
C                                       |      |      |
C                                       | p3 p4| 0  w |
C                                        -------------
C
C
C                to the form, i.e:
C
C
C
C                                        ------
C                                       | w  0 |
C                                       |      |
C                                       | 0  w |
C                   R(T) * P * R  =      -------------
C                                       | a  0 | w  0 |
C                                       |      |      |
C                                       | 0  b | 0  w |
C                                        -------------
C
C
C                in which the offdiagonal block shows now two zeros.
C                Inserting the form of R into the last equation we
C                obtain the matrix equation:
C
C
C                    | a  0 |     |C2 -S2|   | p1 p2|   |C1  S1|
C                    |      |  =  |      | * |      | * |      |
C                    | 0  b |     |S2  C2|   | p3 p4|   |-S1 C1|
C
C
C                from which we obtain two equations in the two unknowns
C                t1 = S1/C1 and t2 = S2/C2:
C
C
C                        t2*p1 - t1*t2*p2 + p3 - t1*p4 = 0
C
C                        t1*p1 - t1*t2*p3 + p2 - t2*p4 = 0
C
C
C                leading to a quadratic equation for both t1 and t2
C                cases:
C
C                              t*t + R*t - 1 = 0
C
C                where:
C
C
C                              p2*p2 + p4*p4 - p1*p1 - p3*p3
C                 R for t1  =  -----------------------------
C                                      p1*p2 + p3*p4
C
C
C                              p3*p3 + p4*p4 - p1*p1 - p2*p2
C                 R for t2  =  -----------------------------
C                                      p1*p3 + p2*p4
C
C
C                Solutions are then given by:
C
C
C                               - R +- sqrt (R*R + 4)
C                          t =  ---------------------
C                                         2
C
C
C                and from the relations t = s/c and s*s + c*c = 1
C                we obtain for the cosines:
C
C
C                       c = sqrt ( 2 / (R*R + 4 -+ R*sqrt (R*R+4) ) )
C
C
C                and, after trying out the solutions, for the sines:
C
C
C                       s = +- sqrt (1 - c*c)
C
C
C                Note in particular the -+ and +- possibilities for
C                the cosines and sines, respectively. Hence a total
C                of 4 sign combinations are possible:
C
C                         For C1,S1   =>   -,+ and +,-
C                         For C2,S2   =>   -,+ and +,-
C
C                corresponding to the 4 different ways the two lobes
C                of the pair (a,b) of two-fold weight degenerate atomic
C                NBOs can be superimposed in space (leave lobes of a
C                fixed and rotate lobes of b by successive 90 degrees):
C
C
C
C
C                                     +(a)+(b)
C
C                                         *
C                                         *
C                           +(a)+(b) *********** -(a)-(b)
C                                         *
C                                         *
C
C                                     -(a)-(b)
C
C
C                 ..................................................
C
C
C                                     +(a)+(b)
C
C                                         *
C                                         *
C                           +(a)-(b) *********** -(a)+(b)
C                                         *
C                                         *
C
C                                     -(a)-(b)
C
C
C                 ..................................................
C
C
C                                     +(a)-(b)
C
C                                         *
C                                         *
C                           +(a)-(b) *********** -(a)+(b)
C                                         *
C                                         *
C
C                                     -(a)+(b)
C
C
C                 ..................................................
C
C
C                                     +(a)-(b)
C
C                                         *
C                                         *
C                           +(a)+(b) *********** -(a)-(b)
C                                         *
C                                         *
C
C                                     -(a)+(b)
C
C
C
C
C                These 4 different cases will lead to the 4 possible
C                (not in correspondence with the above sequence of
C                overlapping lobes pictures) offdiagonal blocks of the
C                occupation matrix P:
C
C
C                    | a  0 |   | b  0 |   | 0  b |   | 0  a |
C                    |      | , |      | , |      | , |      |
C                    | 0  b |   | 0  a |   | a  0 |   | b  0 |
C
C
C                For the sake of uniqueness, we will decide on the
C                lobe overlapping solution that leads to the following
C                offdiagonal block:
C
C
C                            | a  0 |
C                            |      |  and condition  |a| >= |b|
C                            | 0  b |
C
C
C                which has to be found by placing all possible four
C                solutions for C1,S1 and C2,S2 into the equation:
C
C
C                    | a  0 |     |C2 -S2|   | p1 p2|   |C1  S1|
C                    |      |  =  |      | * |      | * |      |
C                    | 0  b |     |S2  C2|   | p3 p4|   |-S1 C1|
C
C
C                and picking the one we are after.
C
C                Relations among the p1,p2,p3,p4 values indicate, if
C                there is the unique set of 4 possible rotations or
C                if there is in fact an infinite possibility of
C                rotations. We have the following two cases:
C
C                  Case i) p2 = -p3 , p1 = p4 or p2 = p3 , p1 = -p4
C                  ------------------------------------------------
C                       This is the case if there is no tilt
C                       angle between the lobe pairs. An infinite
C                       number of rotations is possible here.
C                       Both R values for the t1 and t2 solutions
C                       are given by the indetermined values 0/0.
C                       We proceed here by fixing the first lobe
C                       pairs (i.e. t1 = 0, no angle of rotation)
C                       and rotate the second lobe pair (i.e.
C                       finding the solution for t2 under the
C                       constraint t1 = 0) to match the first lobe
C                       pair.
C
C                  Case ii) p1,p2,p3,p4 not obeing case i)
C                  ---------------------------------------
C                       This is the case if there exists a tilt
C                       angle between the lobe pairs. Proper
C                       solutions exist for both t1 and t2 and
C                       only 4 rotations are possible.
C
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    P            =  full NBAS x NBAS occupation matrix.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    X1,X2        =  NBAS x 2 sections of the NAO
C                                    coefficient matrix in AO basis
C                                    corresponding to the two-fold
C                                    weight degenerate atomic NBOs
C                                    before rotation.
C
C
C                  Output:
C
C                    C1,S1        =  cosine and sine rotation values
C                                    for the first two-fold weight
C                                    degenerate atomic NBO.
C                    C2,S2        =  cosine and sine rotation values
C                                    for the second two-fold weight
C                                    degenerate atomic NBO.
C                    X1,X2        =  NBAS x 2 rotated sections of the
C                                    NAO coefficient matrix in AO basis
C                                    corresponding to the two-fold
C                                    weight degenerate atomic NBOs.
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

         LOGICAL     CASEI
C         LOGICAL     P14ZERO,P23ZERO
         LOGICAL     SAVEC,SAVEP

         INTEGER     N
         INTEGER     NBAS

         DOUBLE PRECISION  A1,A2,A3,A4
         DOUBLE PRECISION  AMAX,BMAX
         DOUBLE PRECISION  C1,C2
         DOUBLE PRECISION  COSA,COSB,COSC,COSD
         DOUBLE PRECISION  ERROR
         DOUBLE PRECISION  F2,F3,F4,F5,F6,F7,F8,F9,F10,F11
         DOUBLE PRECISION  LARGE,VSMALL,EPSILON,EXTREM
         DOUBLE PRECISION  S1,S2
         DOUBLE PRECISION  SINA,SINB,SINC,SIND
         DOUBLE PRECISION  O1,O2,O3,O4
         DOUBLE PRECISION  P1,P2,P3,P4
         DOUBLE PRECISION  R,RNOM,RDENOM
         DOUBLE PRECISION  T2
         DOUBLE PRECISION  Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11
         DOUBLE PRECISION  ZERO,ONE,TWO,FOUR

         DOUBLE PRECISION  XVEC   (1:NBAS)

         DOUBLE PRECISION  X1     (1:NBAS,1:2   )
         DOUBLE PRECISION  X2     (1:NBAS,1:2   )
         DOUBLE PRECISION  P      (1:NBAS,1:NBAS)
         DOUBLE PRECISION  XMAT   (1:NBAS,1:2   )

         PARAMETER   (ZERO    = 0.D0 )
         PARAMETER   (ONE     = 1.D0 )
         PARAMETER   (TWO     = 2.D0 )
         PARAMETER   (FOUR    = 4.D0 )
         PARAMETER   (LARGE   = 1.D+2)
         PARAMETER   (VSMALL  = 1.D-12)
         PARAMETER   (EPSILON = 1.D-12)
         PARAMETER   (EXTREM  = 1.D+12)

         PARAMETER   (F2  =     1.D0 / 2.D0  )
         PARAMETER   (F3  =     3.D0 / 2.D0  )
         PARAMETER   (F4  =    11.D0 / 8.D0  )
         PARAMETER   (F5  =    31.D0 / 8.D0  )
         PARAMETER   (F6  =    69.D0 / 16.D0 )
         PARAMETER   (F7  =   187.D0 / 16.D0 )
         PARAMETER   (F8  =  1843.D0 / 128.D0)
         PARAMETER   (F9  =  4859.D0 / 128.D0)
         PARAMETER   (F10 = 12767.D0 / 256.D0)
         PARAMETER   (F11 = 32965.D0 / 256.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate initial offdiagonal occupation matrix block
C                corresponding to both two-fold weight degenerate
C                atomic NBOs.
C
C
         SAVEC = .TRUE.
         SAVEP = .TRUE.

         CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +              ( NBAS,2,
     +                NBAS,NBAS,
     +                NBAS,2,
     +                NBAS,2,
     +                NBAS,
     +                2,NBAS,2,
     +                0,0,
     +                SAVEC,SAVEP,SAVEC,
     +                X2,P,X1,
     +                XVEC,
     +
     +                           XMAT )
     +
     +
         CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +              ( 6,
     +                ' P sections before ',
     +                NBAS,2,
     +                2,2,
     +                XMAT )
     +
     +
         P1 = XMAT (1,1)
         P2 = XMAT (1,2)
         P3 = XMAT (2,1)
         P4 = XMAT (2,2)
C
C
C             ...check for case i) orientation (lobes only relatively
C                reoriented to each other).
C
C
         CASEI =    (ABS (P2+P3).LT.VSMALL).AND.(ABS (P1-P4).LT.VSMALL)
     +          .OR.(ABS (P2-P3).LT.VSMALL).AND.(ABS (P1+P4).LT.VSMALL)

         IF (CASEI) THEN

             C1 = ONE
             S1 = ZERO

             IF (P4.EQ.ZERO) THEN
                 COSA = ZERO
                 SINA = ONE
             ELSE
                 T2 = P2 / P4
                 COSA = DSQRT (ONE / (ONE + T2*T2))
                 SINA = DSQRT (ONE - COSA * COSA)
             END IF

             COSB = COSA
             COSC = - COSA
             COSD = - COSA
             SINB = - SINA
             SINC = SINA
             SIND = - SINA
C
C
C             ...find largest solution for a in offdiagonal occupation
C                matrix block:
C
C
C                    | a  0 |     |C2 -S2|   | p1 p2|   | 1  0 |
C                    |      |  =  |      | * |      | * |      |
C                    | 0  b |     |S2  C2|   | p3 p4|   | 0  1 |
C
C
C                by forming the 4 possible a1,a2,a3,a4 values using
C                the sines and cosines evaluated above:
C
C
C                               |    C2     S2
C                         -----------------------
C                           A1  |   COSA   SINA
C                           A2  |   COSB   SINB
C                           A3  |   COSC   SINC
C                           A4  |   COSD   SIND
C
C
C                From these 4 possible solutions we exclude those
C                which lead to nonzero offdiagonal elements. The
C                remaining solutions with a positive can either have
C                b = a or b = -a. The latter case is unacceptable
C                and thus an inversion on the second lobe is performed
C                by multiplying the the rotation matrix forming:
C
C
C                        |C2  S2|     | 1  0 |   |C2 -S2|
C                        |      |  =  |      | * |      |
C                        |S2 -C2|     | 0 -1 |   |S2  C2|
C
C
             A1 = COSA * P1 - SINA * P3
             A2 = COSB * P1 - SINB * P3
             A3 = COSC * P1 - SINC * P3
             A4 = COSD * P1 - SIND * P3
             O1 = COSA * P2 - SINA * P4
             O2 = COSB * P2 - SINB * P4
             O3 = COSC * P2 - SINC * P4
             O4 = COSD * P2 - SIND * P4

             IF (ABS (O1).GT.VSMALL) A1 = ZERO
             IF (ABS (O2).GT.VSMALL) A2 = ZERO
             IF (ABS (O3).GT.VSMALL) A3 = ZERO
             IF (ABS (O4).GT.VSMALL) A4 = ZERO

             AMAX = MAX (A1,A2,A3,A4)

             IF (A1.EQ.AMAX) THEN
                 C2 = COSA
                 S2 = SINA
             ELSE IF (A2.EQ.AMAX) THEN
                 C2 = COSB
                 S2 = SINB
             ELSE IF (A3.EQ.AMAX) THEN
                 C2 = COSC
                 S2 = SINC
             ELSE
                 C2 = COSD
                 S2 = SIND
             END IF

             BMAX = S2 * P2 + C2 * P4
C
C
C             ...rotate (possibly with second lobe inversion) second
C                set of degenerate NBOs relative to first set.
C
C
             IF (BMAX.GT.ZERO) THEN

                 IF (ABS (AMAX - BMAX).GT.VSMALL) THEN
                     WRITE (*,*) ' Problems rotating opposite pairs! '
                     WRITE (*,*) ' AMAX,BMAX = ',AMAX,BMAX
                     WRITE (*,*) ' nlo__rotate_pair_degenerate_nbo '
                     WRITE (1,*) ' Problems rotating opposite pairs! '
                     WRITE (1,*) ' AMAX,BMAX = ',AMAX,BMAX
                     WRITE (1,*) ' nlo__rotate_pair_degenerate_nbo '
                     STOP
                 END IF

                 DO N = 1,NBAS
                    R        = C2 * X2 (N,1) - S2 * X2 (N,2)
                    X2 (N,2) = C2 * X2 (N,2) + S2 * X2 (N,1)
                    X2 (N,1) = R
                 END DO

             ELSE

                 IF (ABS (AMAX + BMAX).GT.VSMALL) THEN
                     WRITE (*,*) ' Problems rotating opposite pairs! '
                     WRITE (*,*) ' AMAX,BMAX = ',AMAX,BMAX
                     WRITE (*,*) ' nlo__rotate_pair_degenerate_nbo '
                     WRITE (1,*) ' Problems rotating opposite pairs! '
                     WRITE (1,*) ' AMAX,BMAX = ',AMAX,BMAX
                     WRITE (1,*) ' nlo__rotate_pair_degenerate_nbo '
                     STOP
                 END IF

                 DO N = 1,NBAS
                    R        =   C2 * X2 (N,1) - S2 * X2 (N,2)
                    X2 (N,2) = - C2 * X2 (N,2) - S2 * X2 (N,1)
                    X2 (N,1) = R
                 END DO

             END IF
C
C
C             ...check the newly rotated NBOs.
C
C
             CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                  ( NBAS,2,
     +                    NBAS,NBAS,
     +                    NBAS,2,
     +                    NBAS,2,
     +                    NBAS,
     +                    2,NBAS,2,
     +                    0,0,
     +                    SAVEC,SAVEP,SAVEC,
     +                    X2,P,X1,
     +                    XVEC,
     +
     +                               XMAT )
     +
     +
            CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                 ( 6,
     +                   ' P sections after relative reorientation ',
     +                   NBAS,2,
     +                   2,2,
     +                   XMAT )
     +
     +
            RETURN

         END IF
C
C
C             ...the absolute reorientation case ii).
C                Find the two sine and cosine solutions for rotation
C                of first degenerate NBOs. If either pair P1,P4 or
C                P2,P3 values are below the calculation accuracy of
C                VSMALL, we apply no rotation. This catch is necessary
C                as the other values might be very small as well,
C                leading to a false conclusion of necessity of small
C                rotations.
C
C
C         P14ZERO = ABS (P1).LT.VSMALL .AND. ABS (P4).LT.VSMALL
C         P23ZERO = ABS (P2).LT.VSMALL .AND. ABS (P3).LT.VSMALL
C
C         IF (P14ZERO) THEN
C             R = EXTREM * DSIGN (ONE,P2-P3)
C         ELSE IF (P23ZERO) THEN
C             R = EXTREM * DSIGN (ONE,P4-P1)
C         ELSE
C             R = (P2*P2 + P4*P4 - P1*P1 - P3*P3) / (P1*P2 + P3*P4)
C         END IF

         ERROR  = (ABS (P1) + ABS (P2) + ABS (P3) + ABS (P4)) * EPSILON
         RNOM   = P2*P2 + P4*P4 - P1*P1 - P3*P3
         RDENOM = P1*P2 + P3*P4

         IF (ABS (RDENOM).LT.ERROR) THEN
             R = EXTREM * DSIGN (ONE,RNOM)
         ELSE
             R = RNOM / RDENOM
         END IF

         WRITE (*,*) ' R for A,B = ',R

         IF (ABS (R).GE.EXTREM) THEN

             IF (R.LT.ZERO) THEN
                 COSA = ONE
                 SINB = ONE
                 SINA = ZERO
                 COSB = ZERO
             ELSE
                 COSB = ONE
                 SINA = - ONE
                 SINB = ZERO
                 COSA = ZERO
             END IF

         ELSE IF (ABS (R).GT.LARGE) THEN

             Z1  = ONE / ABS (R)
             Z2  = Z1  * Z1
             Z3  = Z2  * Z1
             Z4  = Z3  * Z1
             Z5  = Z4  * Z1
             Z6  = Z5  * Z1
             Z7  = Z6  * Z1
             Z8  = Z7  * Z1
             Z9  = Z8  * Z1
             Z10 = Z9  * Z1
             Z11 = Z10 * Z1

             IF (R.LT.ZERO) THEN
                 COSA = ONE - F2*Z2 + F4*Z4 - F6*Z6 + F8*Z8 - F10*Z10
                 COSB =  Z1 - F3*Z3 + F5*Z5 - F7*Z7 + F9*Z9 - F11*Z11
                 SINA = - COSB
                 SINB = COSA
             ELSE
                 COSB = ONE - F2*Z2 + F4*Z4 - F6*Z6 + F8*Z8 - F10*Z10
                 SINB =  Z1 - F3*Z3 + F5*Z5 - F7*Z7 + F9*Z9 - F11*Z11
                 SINA = - COSB
                 COSA = SINB
             END IF

         ELSE

             Z1 = R * R + FOUR
             Z2 = R * DSQRT (Z1)
             Z3 = TWO / (Z1 + Z2)
             Z4 = TWO / (Z1 - Z2)
             COSA = DSQRT (Z3)
             COSB = DSQRT (Z4)
             SINA = - DSQRT (ONE - Z3)
             SINB =   DSQRT (ONE - Z4)

         END IF
C
C
C             ...find the two sine and cosine solutions for rotation
C                of second degenerate NBOs.
C
C
C         IF (P14ZERO) THEN
C             R = EXTREM * DSIGN (ONE,P3-P2)
C         ELSE IF (P23ZERO) THEN
C             R = EXTREM * DSIGN (ONE,P4-P1)
C         ELSE
C             R = (P3*P3 + P4*P4 - P1*P1 - P2*P2) / (P1*P3 + P2*P4)
C         END IF

         ERROR  = (ABS (P1) + ABS (P2) + ABS (P3) + ABS (P4)) * EPSILON
         RNOM   = P3*P3 + P4*P4 - P1*P1 - P2*P2
         RDENOM = P1*P3 + P2*P4

         IF (ABS (RDENOM).LT.ERROR) THEN
             R = EXTREM * DSIGN (ONE,RNOM)
         ELSE
             R = RNOM / RDENOM
         END IF

         WRITE (*,*) ' R for C,D = ',R

         IF (ABS (R).GE.EXTREM) THEN

             IF (R.LT.ZERO) THEN
                 COSC = ONE
                 COSD = ZERO
                 SINC = ZERO
                 SIND = ONE
             ELSE
                 COSD = ONE
                 SIND = ZERO
                 SINC = - ONE
                 COSC = ZERO
             END IF

         ELSE IF (ABS (R).GT.LARGE) THEN

             Z1  = ONE / ABS (R)
             Z2  = Z1  * Z1
             Z3  = Z2  * Z1
             Z4  = Z3  * Z1
             Z5  = Z4  * Z1
             Z6  = Z5  * Z1
             Z7  = Z6  * Z1
             Z8  = Z7  * Z1
             Z9  = Z8  * Z1
             Z10 = Z9  * Z1
             Z11 = Z10 * Z1

             IF (R.LT.ZERO) THEN
                 COSC = ONE - F2*Z2 + F4*Z4 - F6*Z6 + F8*Z8 - F10*Z10
                 COSD =  Z1 - F3*Z3 + F5*Z5 - F7*Z7 + F9*Z9 - F11*Z11
                 SINC = - COSD
                 SIND = COSC
             ELSE
                 COSD = ONE - F2*Z2 + F4*Z4 - F6*Z6 + F8*Z8 - F10*Z10
                 SIND =  Z1 - F3*Z3 + F5*Z5 - F7*Z7 + F9*Z9 - F11*Z11
                 SINC = - COSD
                 COSC = SIND
             END IF

         ELSE

             Z1 = R * R + FOUR
             Z2 = R * DSQRT (Z1)
             Z3 = TWO / (Z1 + Z2)
             Z4 = TWO / (Z1 - Z2)
             COSC = DSQRT (Z3)
             COSD = DSQRT (Z4)
             SINC = - DSQRT (ONE - Z3)
             SIND =   DSQRT (ONE - Z4)

         END IF

         WRITE (*,*) ' COSA = ',COSA
         WRITE (*,*) ' SINA = ',SINA
         WRITE (*,*) ' COSB = ',COSB
         WRITE (*,*) ' SINB = ',SINB
         WRITE (*,*) ' COSC = ',COSC
         WRITE (*,*) ' SINC = ',SINC
         WRITE (*,*) ' COSD = ',COSD
         WRITE (*,*) ' SIND = ',SIND

         WRITE (*,*) ' A ident = ',COSA*COSA + SINA*SINA
         WRITE (*,*) ' B ident = ',COSB*COSB + SINB*SINB
         WRITE (*,*) ' C ident = ',COSC*COSC + SINC*SINC
         WRITE (*,*) ' D ident = ',COSD*COSD + SIND*SIND
C
C
C             ...find largest solution for |a| in offdiagonal occupation
C                matrix block:
C
C
C                    | a  0 |     |C2 -S2|   | p1 p2|   |C1  S1|
C                    |      |  =  |      | * |      | * |      |
C                    | 0  b |     |S2  C2|   | p3 p4|   |-S1 C1|
C
C
C                by forming the 4 possible a1,a2,a3,a4 values using
C                the sines and cosines evaluated above:
C
C
C                         |    C1     S1     C2     S2
C                         --------------------------------
C                     A1  |   COSA   SINA   COSC   SINC
C                     A2  |   COSA   SINA   COSD   SIND
C                     A3  |   COSB   SINB   COSC   SINC
C                     A4  |   COSB   SINB   COSD   SIND
C
C
         A1 = COSA*COSC*P1 - SINA*COSC*P2 - COSA*SINC*P3 + SINA*SINC*P4
         A2 = COSA*COSD*P1 - SINA*COSD*P2 - COSA*SIND*P3 + SINA*SIND*P4
         A3 = COSB*COSC*P1 - SINB*COSC*P2 - COSB*SINC*P3 + SINB*SINC*P4
         A4 = COSB*COSD*P1 - SINB*COSD*P2 - COSB*SIND*P3 + SINB*SIND*P4

         A1 = DABS (A1)
         A2 = DABS (A2)
         A3 = DABS (A3)
         A4 = DABS (A4)

         AMAX = MAX (A1,A2,A3,A4)

         IF (A1.EQ.AMAX) THEN
             C1 = COSA
             C2 = COSC
             S1 = SINA
             S2 = SINC
             WRITE (*,*) ' Took A,C ! '
         ELSE IF (A2.EQ.AMAX) THEN
             C1 = COSA
             C2 = COSD
             S1 = SINA
             S2 = SIND
             WRITE (*,*) ' Took A,D ! '
         ELSE IF (A3.EQ.AMAX) THEN
             C1 = COSB
             C2 = COSC
             S1 = SINB
             S2 = SINC
             WRITE (*,*) ' Took B,C ! '
         ELSE
             C1 = COSB
             C2 = COSD
             S1 = SINB
             S2 = SIND
             WRITE (*,*) ' Took B,D ! '
         END IF

         AMAX = C2 * (C1*P1 - S1*P2) - S2 * (C1*P3 - S1*P4)
         BMAX = S2 * (S1*P1 + C1*P2) + C2 * (S1*P3 + C1*P4)
C
C
C             ...rotate both sets of degenerate NBOs.
C
C
         IF (AMAX.GT.ZERO) THEN
             DO N = 1,NBAS
                R        = C1 * X1 (N,1) - S1 * X1 (N,2)
                X1 (N,2) = C1 * X1 (N,2) + S1 * X1 (N,1)
                X1 (N,1) = R
             END DO
         ELSE
             DO N = 1,NBAS
                R        = - C1 * X1 (N,1) + S1 * X1 (N,2)
                X1 (N,2) =   C1 * X1 (N,2) + S1 * X1 (N,1)
                X1 (N,1) = R
             END DO
         END IF

         IF (BMAX.GT.ZERO) THEN
             DO N = 1,NBAS
                R        = C2 * X2 (N,1) - S2 * X2 (N,2)
                X2 (N,2) = C2 * X2 (N,2) + S2 * X2 (N,1)
                X2 (N,1) = R
             END DO
         ELSE
             DO N = 1,NBAS
                R        =   C2 * X2 (N,1) - S2 * X2 (N,2)
                X2 (N,2) = - C2 * X2 (N,2) - S2 * X2 (N,1)
                X2 (N,1) = R
             END DO
         END IF
C
C
C             ...check the newly rotated NBOs.
C
C
         CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +              ( NBAS,2,
     +                NBAS,NBAS,
     +                NBAS,2,
     +                NBAS,2,
     +                NBAS,
     +                2,NBAS,2,
     +                0,0,
     +                SAVEC,SAVEP,SAVEC,
     +                X2,P,X1,
     +                XVEC,
     +
     +                           XMAT )
     +
     +
         CALL  MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +              ( 6,
     +                ' P sections after ',
     +                NBAS,2,
     +                2,2,
     +                XMAT )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
