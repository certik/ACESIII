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
         SUBROUTINE  NLO__FORM_SP2_PAIR_NBO
     +
     +                    ( NBAS,
     +                      P,
     +                      XVEC,
     +
     +                            HA,HB,
     +                            SA,PA,
     +                            SB,PB )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_SP2_PAIR_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms a pair of linear sp2 combinations
C                from two sets of three input NBOs (1 totally symmetric
C                s-type NBO and two p-type NBOs). The sp2 combinations
C                formed are as follows:
C
C
C                         ha1  =  m * sa + 2n * pax
C                         ha2  =  m * sa -  n * pax + q * pay
C                         ha3  =  m * sa -  n * pax - q * pay
C
C                         hb1  =  m * sb + 2n * pbx
C                         hb2  =  m * sb -  n * pbx + q * pby
C                         hb3  =  m * sb -  n * pbx - q * pby
C
C                where:
C                                 m = 1 / sqrt (3)
C                                 n = 1 / sqrt (6)
C                                 q = 1 / sqrt (2)
C
C
C                The final pair of sp2 combinations is such that their
C                main directional lobes corresponding to ha1 and hb1
C                are either eclipsed or in anti position. How is that
C                achieved? Suppose the incomming px and py functions
C                of both input NBO sets are already ordered such that
C                their lobes match and form rectangular angles between
C                them (px is labeled as 'x' and py is labeled as '*'):
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
C                For each of the px,py pairs there are thus 4 different
C                ways to form the main directional lobe hybrid (ha1 and
C                hb1), depending on where the sp2 lobe will point to:
C                +x axis, -x axis, +y axis and -y axis. Mathematically
C                this would correspond to a sign change and px <-> py
C                interchange when forming ha1 (or hb1):
C
C
C                     ha1  =  m * sa + 2n * pax
C                     ha1  =  m * sa - 2n * pax   (sign)
C                     ha1  =  m * sa + 2n * pyx   (px <-> py)
C                     ha1  =  m * sa - 2n * pyx   (sign and px <-> py)
C
C
C                We now form the occupation matrix element PH between
C                both ha1 and hb1. Of the total of 16 possible values
C                (4 x 4 for the 4 ha1 and 4 hb1 orientations), only
C                those corresponding to the eclipsed and anti positions
C                between ha1 and hb1 will occur only once. To get a
C                unique orientation we will thus look for the largest
C                absolute non-degenerate value among all 16 possible
C                PH values. Once identified we form the rest of the
C                corresponding sp2 hybrids (ha2,ha3,hb2,hb3).
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    P            =  full NBAS x NBAS occupation matrix
C                                    in AO basis.
C                    XVEC         =  flp scratch array of vector type.
C                    SA           =  s-type NBO coefficient vector in
C                                    AO basis for first NBO triplet.
C                    PA           =  NBAS x 2 px- and py-type NBO
C                                    coefficient matrix in AO basis
C                                    for first NBO triplet.
C                    SB           =  s-type NBO coefficient vector in
C                                    AO basis for second NBO triplet.
C                    PB           =  NBAS x 2 px- and py-type NBO
C                                    coefficient matrix in AO basis
C                                    for second NBO triplet.
C
C
C                  Output:
C
C                    HA           =  3 x 3 matrix, containing the
C                                    sp2 hybridization coefficients:
C                                    1st column defines ha1, 2nd column
C                                    defines ha2, 3rd column defines ha3
C                                    in basis (rows): sa,pax,pay
C                    HB           =  3 x 3 matrix, containing the
C                                    sp2 hybridization coefficients:
C                                    1st column defines hb1, 2nd column
C                                    defines hb2, 3rd column defines hb3
C                                    in basis (rows): sb,pbx,pby
C                    SA           =  main directional lobe hybrid NBO
C                                    coefficient vector in AO basis
C                                    for first NBO triplet.
C                    PA           =  NBAS x 2 second and third hybrid
C                                    NBO coefficient matrix in AO basis
C                                    for first NBO triplet.
C                    SB           =  main directional lobe hybrid NBO
C                                    coefficient vector in AO basis
C                                    for second NBO triplet.
C                    PB           =  NBAS x 2 second and third hybrid
C                                    NBO coefficient matrix in AO basis
C                                    for second NBO triplet.
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

         LOGICAL     ABORT
         LOGICAL     ALERT
         LOGICAL     EQUAL

         INTEGER     CASE
         INTEGER     K,L,M,N
         INTEGER     NBAS
         INTEGER     PXA,PXB,PYA,PYB

         INTEGER     HYB  (1:2 )
         INTEGER     POSA (1:16)
         INTEGER     POSB (1:16)

         DOUBLE PRECISION  E,F,G
         DOUBLE PRECISION  MA,NA,QA,MB,NB,QB
         DOUBLE PRECISION  PHK,PHL
         DOUBLE PRECISION  SGA,SGB
         DOUBLE PRECISION  SQR3TH,SQR6TH,SQR2TH
         DOUBLE PRECISION  SYMACC
         DOUBLE PRECISION  X,Y,Z
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  PARA (1:16  )
         DOUBLE PRECISION  PARB (1:16  )
         DOUBLE PRECISION  PH   (1:16  )
         DOUBLE PRECISION  SA   (1:NBAS)
         DOUBLE PRECISION  SB   (1:NBAS)
         DOUBLE PRECISION  XVEC (1:NBAS)

         DOUBLE PRECISION  HA   (1:3,1:3)
         DOUBLE PRECISION  HB   (1:3,1:3)

         DOUBLE PRECISION  P    (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PA   (1:NBAS,1:2   )
         DOUBLE PRECISION  PB   (1:NBAS,1:2   )

         PARAMETER   (ZERO   = 0.D0                                   )
         PARAMETER   (SYMACC = 1.D-10                                 )
         PARAMETER   (SQR3TH = 0.57735026918962576450914878050195750D0)
         PARAMETER   (SQR6TH = 0.40824829046386301636621401245098191D0)
         PARAMETER   (SQR2TH = 0.70710678118654752440084436210484909D0)

         DATA  PARA  /+1.D0,+1.D0,+1.D0,+1.D0,+1.D0,+1.D0,+1.D0,+1.D0,
     +                -1.D0,-1.D0,-1.D0,-1.D0,-1.D0,-1.D0,-1.D0,-1.D0/
         DATA  POSA  /  1  ,  1  ,  1  ,  1  ,  2  ,  2  ,  2  ,  2  ,
     +                  1  ,  1  ,  1  ,  1  ,  2  ,  2  ,  2  ,  2  /
         DATA  PARB  /+1.D0,+1.D0,-1.D0,-1.D0,+1.D0,+1.D0,-1.D0,-1.D0,
     +                +1.D0,+1.D0,-1.D0,-1.D0,+1.D0,+1.D0,-1.D0,-1.D0/
         DATA  POSB  /  1  ,  2  ,  1  ,  2  ,  1  ,  2  ,  1  ,  2  ,
     +                  1  ,  2  ,  1  ,  2  ,  1  ,  2  ,  1  ,  2  /
C
C
C------------------------------------------------------------------------
C
C
C             ...form the 16 occupation matrix elements between
C                the main directional lobe hybrids.
C
C
         DO CASE = 1,16

            SGA = PARA (CASE)
            SGB = PARB (CASE)
            PXA = POSA (CASE)
            PXB = POSB (CASE)

            IF (MOD (CASE,4).EQ.1) THEN
                DO N = 1,NBAS
                   E = SQR3TH * SA (N)
                   F = SGA * SQR6TH * PA (N,PXA)
                   XVEC (N) = E + F + F
                END DO
            END IF

            Z = ZERO
            DO L = 1,NBAS
               X = ZERO
               DO K = 1,NBAS
                  X = X + XVEC (K) * P (K,L)
               END DO
               E = SQR3TH * SB (L)
               F = SGB * SQR6TH * PB (L,PXB)
               Y = E + F + F
               Z = Z + X * Y
            END DO

            PH (CASE) = Z

         END DO
C
C
C             ...identify those non-degenerate occupation matrix
C                elements corresponding to the eclipsed and anti
C                positions of the main directional lobe hybrids.
C
C
         M = 0

         DO K = 1,15
            PHK = PH (K)
            IF (PHK.NE.ZERO) THEN
                EQUAL = .FALSE.
                DO L = K+1,16
                   PHL = PH (L)
                   IF (ABS (PHK-PHL).LT.SYMACC) THEN
                       EQUAL = .TRUE.
                       PH (L) = ZERO
                   END IF
                END DO
                IF (.NOT.EQUAL) THEN
                    M = M + 1
                    HYB (M) = K
                END IF
            END IF
         END DO

         IF (PH (16).NE.ZERO) THEN
             M = M + 1
             HYB (M) = 16
         END IF
C
C
C             ...ideally only 2 out of the 16 occupation matrix
C                elements were found. If none were found, we have
C                an accidental degeneracy problem and we must stop.
C                If more than 2 were found, we have an accuracy
C                problem, but the program can still try to continue.
C
C
         ABORT = M .EQ. 0
         ALERT = M .GT. 2

         IF (ABORT) THEN
             WRITE (*,*) ' Problem in sp2 hybrid orientation! '
             WRITE (*,*) ' Cannot find unique orientation! '
             WRITE (*,*) ' nlo__form_sp2_pair_nbo '
             WRITE (1,*) ' Problem in sp2 hybrid orientation! '
             WRITE (1,*) ' Cannot find unique orientation! '
             WRITE (1,*) ' nlo__form_sp2_pair_nbo '
             STOP
         END IF

         IF (ALERT) THEN
             WRITE (*,*) ' Floppy sp2 hybrid orientation! '
             WRITE (*,*) ' Will proceed ... '
             WRITE (1,*) ' Floppy sp2 hybrid orientation! '
             WRITE (1,*) ' Will proceed ... '
         END IF

C         WRITE (*,*) ' Hybrids = ',(HYB(K),K=1,M)
C
C
C             ...the maximum of the surviving 2 occupation matrix
C                elements is chosen and the corresponding sign
C                changes and/or px <-> py interchanges determined.
C
C
         K = HYB (1)
         L = HYB (2)
         PHK = ABS (PH (K))
         PHL = ABS (PH (L))

         X = MAX (PHK,PHL)
         IF (X.EQ.PHK) THEN
             CASE = K
         ELSE
             CASE = L
         END IF

         SGA = PARA (CASE)
         SGB = PARB (CASE)
         PXA = POSA (CASE)
         PXB = POSB (CASE)

         MA = SQR3TH
         NA = SGA * SQR6TH
         QA = SGA * SQR2TH
         MB = SQR3TH
         NB = SGB * SQR6TH
         QB = SGB * SQR2TH
C
C
C             ...form the hybridization matrices.
C
C
         IF (PXA.EQ.1) THEN
             PYA = 2
             HA (1,1) = MA
             HA (2,1) = NA + NA
             HA (3,1) = ZERO
             HA (1,2) = MA
             HA (2,2) = - NA
             HA (3,2) = QA
             HA (1,3) = MA
             HA (2,3) = - NA
             HA (3,3) = - QA
         ELSE
             PYA = 1
             HA (1,1) = MA
             HA (2,1) = ZERO
             HA (3,1) = NA + NA
             HA (1,2) = MA
             HA (2,2) = QA
             HA (3,2) = - NA
             HA (1,3) = MA
             HA (2,3) = - QA
             HA (3,3) = - NA
         END IF

         IF (PXB.EQ.1) THEN
             PYB = 2
             HB (1,1) = MB
             HB (2,1) = NB + NB
             HB (3,1) = ZERO
             HB (1,2) = MB
             HB (2,2) = - NB
             HB (3,2) = QB
             HB (1,3) = MB
             HB (2,3) = - NB
             HB (3,3) = - QB
         ELSE
             PYB = 1
             HB (1,1) = MB
             HB (2,1) = ZERO
             HB (3,1) = NB + NB
             HB (1,2) = MB
             HB (2,2) = QB
             HB (3,2) = - NB
             HB (1,3) = MB
             HB (2,3) = - QB
             HB (3,3) = - NB
         END IF
C
C
C             ...form the final sp2 hybrid functions.
C
C
         DO N = 1,NBAS
            E = MA * SA (N)
            F = NA * PA (N,PXA)
            G = QA * PA (N,PYA)
            SA (N)   = E + F + F
            PA (N,1) = E - F + G
            PA (N,2) = E - F - G
         END DO

         DO N = 1,NBAS
            E = MB * SB (N)
            F = NB * PB (N,PXB)
            G = QB * PB (N,PYB)
            SB (N)   = E + F + F
            PB (N,1) = E - F + G
            PB (N,2) = E - F - G
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
