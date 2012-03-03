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
         SUBROUTINE  NLO__ANALYZE_RINGS
     +
     +                    ( N,
     +                      R,
     +                      X,Y,Z,
     +                      THRESH,
     +                      USED,
     +
     +                              PLANAR )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ANALYZE_RINGS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine checks, if the current ring R consists
C                of subrings of lesser size.
C
C                Procedure:
C
C                  a) Take the first point P1
C
C
C
C                  Input:
C
C                    N            =  size of ring
C                    R            =  ring atomic indices
C                    X,Y,Z        =  their x,y,z coordinates
C                    THRESH       =  threshold for determining
C                                    coplanarity
C
C
C                  Output:
C
C                    PLANAR       =  is true, if all points are
C                                    coplanar, false otherwise
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

         LOGICAL     LARGE
         LOGICAL     PLANAR

         INTEGER     I,K,N

         INTEGER     R    (1:N)
         INTEGER     USED (1:N)

         DOUBLE PRECISION  A,B,C,D
         DOUBLE PRECISION  DIFF
         DOUBLE PRECISION  M1,MI
         DOUBLE PRECISION  THRESH
         DOUBLE PRECISION  X1,Y1,Z1
         DOUBLE PRECISION  XC,YC,ZC
         DOUBLE PRECISION  XN
         DOUBLE PRECISION  XP,YP,ZP
         DOUBLE PRECISION  XI,YI,ZI
         DOUBLE PRECISION  ZERO,TWOPI

         DOUBLE PRECISION  COSINE,THETA

         DOUBLE PRECISION  F (1:100)
         DOUBLE PRECISION  X (1:N)
         DOUBLE PRECISION  Y (1:N)
         DOUBLE PRECISION  Z (1:N)

         PARAMETER  (TWOPI = 360.D0)
         PARAMETER  (ZERO  = 0.D0  )
#ifdef __GNUC__
         DOUBLE PRECISION DACOSD
         EXTERNAL DACOSD
#endif
C
C
C------------------------------------------------------------------------
C
C
C             ...check for coplanarity of ring members.
C
C
         CALL  NLO__CHECK_COPLANARITY
     +
     +              ( N,
     +                X,Y,Z,
     +                THRESH,
     +
     +                        PLANAR,
     +                        A,B,C,D )
     +
     +
C
C
C             ...if coplanar, decompose points into regular polygons.
C
C
         IF (PLANAR) THEN

             XC = ZERO
             YC = ZERO
             ZC = ZERO
             DO I = 1,N
                XC = XC + X (I)
                YC = YC + Y (I)
                ZC = ZC + Z (I)
                USED (I) = 0
             END DO
             XC = XC / DFLOAT (N)
             YC = YC / DFLOAT (N)
             ZC = ZC / DFLOAT (N)

             XN = DFLOAT (N)

             X1 = X (1) - XC
             Y1 = Y (1) - YC
             Z1 = Z (1) - ZC
             M1 = DSQRT (X1 * X1 + Y1 * Y1 + Z1 * Z1)
             F (1) = ZERO

             DO I = 2,N
                XI = X (I) - XC
                YI = Y (I) - YC
                ZI = Z (I) - ZC
                MI = DSQRT (XI * XI + YI * YI + ZI * ZI)
                XP = Y1 * ZI - Z1 * YI
                YP = Z1 * XI - X1 * ZI
                ZP = X1 * YI - Y1 * XI
                LARGE = (A*XP + B*YP + C*ZP) .LT. ZERO
                COSINE = (X1*XI + Y1*YI + Z1*ZI) / (M1*MI)
                THETA = DACOSD (COSINE)
                IF (LARGE) THEN
                    THETA = TWOPI - THETA
                END IF
                F (I) = XN * THETA / TWOPI
                K = NINT (F (I))
                DIFF = ABS (F (I) - DFLOAT (K))
                WRITE (*,*) ' Theta,F,DIFF = ',THETA,F(I),DIFF
             END DO



         END IF




C             DO I = 1,1
C                IF (USED (I).EQ.0) THEN
C                    XI = X (I) - XC
C                    YI = Y (I) - YC
C                    ZI = Z (I) - ZC
C                    MI = DSQRT (XI * XI + YI * YI + ZI * ZI)
C
C                    NSUB = 1
C                    RSUB (1) = R (I)
C
C                    DO J = I+1,N
C                       IF (USED (J).EQ.0) THEN
C                           XJ = X (J) - XC
C                           YJ = Y (J) - YC
C                           ZJ = Z (J) - ZC
C                           MJ = DSQRT (XJ * XJ + YJ * YJ + ZJ * ZJ)
C                           XP = YI * ZJ - ZI * YJ
C                           YP = ZI * XJ - XI * ZJ
C                           ZP = XI * YJ - YI * XJ
C                           LARGE = (A*XP + B*YP + C*ZP) .LT. ZERO
C                           COSINE = (XI*XJ + YI*YJ + ZI*ZJ) / (MI*MJ)
C                           THETA = DACOSD (COSINE)
C                           IF (LARGE) THEN
C                               THETA = TWOPI - THETA
C                           END IF
C                           F = XN * THETA / TWOPI
C                           K = NINT (F)
C                           DIFF = ABS (F - DFLOAT (K))
C                           WRITE (*,*) ' Theta,F,DIFF = ',THETA,F,DIFF
C                       END IF
C                    END DO
C                END IF
C             END DO
C         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
