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
         SUBROUTINE  NLO__FORM_SP2_NBO
     +
     +                    ( NBAS,
     +
     +                            H11,H12,H13,
     +                            H21,H22,H23,
     +                            H31,H32,H33,
     +                            S,P )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_SP2_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms trigonal sp2 hybrids from one
C                totally symmetric s-type NBO and two p-type NBOs.
C                The sp2 hybrids are formed as follows:
C
C
C                         h1  =  m * s + 2n * p1
C                         h2  =  m * s -  n * p1 + q * p2
C                         h3  =  m * s -  n * p1 - q * p2
C
C                where:
C                                 m = 1 / sqrt (3)
C                                 n = 1 / sqrt (6)
C                                 q = 1 / sqrt (2)
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NHYB         =  total # of AO's in AO basis
C                    S            =  s-type NBO coefficient vector in
C                                    AO basis.
C                    P            =  NBAS x 2 p1- and p2-type NBO
C                                    coefficient matrix in AO basis.
C
C
C                  Output:
C
C                    Hij          =  3 x 3 matrix ij components,
C                                    containing the sp2 hybridization
C                                    coefficients: elements H11,H21,H31
C                                    define h1, elements H12,H22,H32
C                                    define h2, elements H13,H23,H33
C                                    define h3 in basis (rows): s,p1,p2
C                    S            =  h1 hybrid NBO coefficient vector
C                                    in AO basis.
C                    P            =  h2 and h3 hybrid NBO coefficient
C                                    vectors in AO basis.
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

         INTEGER     N
         INTEGER     NBAS

         DOUBLE PRECISION  A,B,C
         DOUBLE PRECISION  H11,H12,H13,H21,H22,H23,H31,H32,H33
         DOUBLE PRECISION  SQR3TH,SQR6TH,SQR2TH
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  S  (1:NBAS)
         DOUBLE PRECISION  P  (1:NBAS,1:2)

         PARAMETER   (ZERO   = 0.D0)
         PARAMETER   (SQR3TH = 0.57735026918962576450914878050195750D0)
         PARAMETER   (SQR6TH = 0.40824829046386301636621401245098191D0)
         PARAMETER   (SQR2TH = 0.70710678118654752440084436210484909D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...form the sp2 hybridization matrix elements.
C
C
         H11 = SQR3TH
         H21 = SQR6TH + SQR6TH
         H31 = ZERO
         H12 = SQR3TH
         H22 = - SQR6TH
         H32 = SQR2TH
         H13 = SQR3TH
         H23 = - SQR6TH
         H33 = - SQR2TH
C
C
C             ...form the sp2 hybrid functions.
C
C
         DO N = 1,NBAS
            A = SQR3TH * S (N)
            B = SQR6TH * P (N,1)
            C = SQR2TH * P (N,2)
            S (N)   = A + B + B
            P (N,1) = A - B + C
            P (N,2) = A - B - C
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
