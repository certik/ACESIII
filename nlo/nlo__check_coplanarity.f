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
         SUBROUTINE  NLO__CHECK_COPLANARITY
     +
     +                    ( N,
     +                      X,Y,Z,
     +                      THRESH,
     +
     +                              PLANAR,
     +                              A,B,C,D )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__CHECK_COPLANARITY
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine checks, if a set of N points is coplanar
C                in space. The location of the points are given by
C                their x,y,z coordinates. If the points are found to be
C                coplanar the routine also returns the coefficients
C                a,b,c,d that define the plane ax+by+cz+d=0, where
C                u=(a,b,c) is the normal vector of the plane with a,b,c
C                being the normalized components: a*a+b*b+c*c=1.
C
C                Procedure:
C                
C                  a) Take the first point P1
C
C                  b) Find a pair of points PI,PJ from the set of N-1
C                     remaining points, such that the resulting cross
C                     product vector CPV between vectors P1->PI and
C                     P1->PJ has largest magnitude. If largest magnitude
C                     is lower than the zero threshold value all points
C                     are considered to be colinear and thus also
C                     coplanar.
C
C                  c) Calculate the scalar products SPR between the CPV
C                     found in b) and all remaining vectors P1->PK for
C                     the N-3 remaining K points. Form the maximum
C                     SPR (max) of all the SPRs.
C
C                  d) If SPR (max) < zero threshold, all points are
C                     considered to be coplanar.
C
C
C                  Input:
C
C                    N            =  # of points to be checked
C                    X,Y,Z        =  their x,y,z coordinates
C                    THRESH       =  zero threshold, determining
C                                    coplanarity
C
C
C                  Output:
C
C                    PLANAR       =  is true, if all points are
C                                    coplanar, false otherwise
C                    A,B,C,D      =  if planar is true, these are
C                                    the coefficients defining the
C                                    plane where all points ly on.
C                                    If planar is false, they are
C                                    set equal to zero.
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

         LOGICAL     PLANAR

         INTEGER     I,J,K,L,N

         DOUBLE PRECISION  A,B,C,D
         DOUBLE PRECISION  MAG,MAXMAG
         DOUBLE PRECISION  SPR,MAXSPR
         DOUBLE PRECISION  THRESH
         DOUBLE PRECISION  X1,Y1,Z1
         DOUBLE PRECISION  XCP,YCP,ZCP
         DOUBLE PRECISION  XVI,YVI,ZVI,XVJ,YVJ,ZVJ
         DOUBLE PRECISION  XVK,YVK,ZVK,XVL,YVL,ZVL
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  X (1:N)
         DOUBLE PRECISION  Y (1:N)
         DOUBLE PRECISION  Z (1:N)

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...deal with trivial cases.
C
C
         A = ZERO
         B = ZERO
         C = ZERO
         D = ZERO

         IF (N.LT.3) THEN
             PLANAR = .FALSE.
             RETURN
         END IF
C
C
C             ...find vectors PI and PJ.
C
C
         X1 = X (1)
         Y1 = Y (1)
         Z1 = Z (1)

         MAXMAG = ZERO

         DO K = 2,N
            XVK = X (K) - X1
            YVK = Y (K) - Y1
            ZVK = Z (K) - Z1
            DO L = K+1,N
               XVL = X (L) - X1
               YVL = Y (L) - Y1
               ZVL = Z (L) - Z1
               XCP = YVK * ZVL - ZVK * YVL
               YCP = ZVK * XVL - XVK * ZVL
               ZCP = XVK * YVL - YVK * XVL
               MAG = DSQRT (XCP * XCP + YCP * YCP + ZCP * ZCP)
               IF (MAG.GT.MAXMAG) THEN
                   I = K
                   J = L
                   MAXMAG = MAG
               END IF
            END DO
         END DO
C
C
C             ...form cross product components between P1->PI and
C                P1->PJ.
C
C
         XVI = X (I) - X1
         YVI = Y (I) - Y1
         ZVI = Z (I) - Z1
         XVJ = X (J) - X1
         YVJ = Y (J) - Y1
         ZVJ = Z (J) - Z1
         XCP = YVI * ZVJ - ZVI * YVJ
         YCP = ZVI * XVJ - XVI * ZVJ
         ZCP = XVI * YVJ - YVI * XVJ
C
C
C             ...form maximum of scalar products between the cross
C                product vector and the remaining N-3 vectors P1->PK.
C
C
         MAXSPR = ZERO

         DO K = 2,N
            IF ((K.NE.I).AND.(K.NE.J)) THEN
                 XVK = X (K) - X1
                 YVK = Y (K) - Y1
                 ZVK = Z (K) - Z1
                 SPR = XCP * XVK + YCP * YVK + ZCP * ZVK
                 MAXSPR = MAX (ABS (SPR),MAXSPR)
            END IF
         END DO
C
C
C             ...decide if set of points is coplanar and if they
C                are calculate the coefficients a,b,c,d that define
C                the equation of the plane ax+by+cz+d=0, where
C                u=(a,b,c) defines the normal vector of the plane.
C
C
         PLANAR = MAXSPR .LT. THRESH

         IF (PLANAR) THEN
             MAG = DSQRT (XCP * XCP + YCP * YCP + ZCP * ZCP)
             A = XCP / MAG
             B = YCP / MAG
             C = ZCP / MAG
             D = - (XCP * X1 + YCP * Y1 + ZCP * Z1) / MAG
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
