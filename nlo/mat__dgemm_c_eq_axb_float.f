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
         SUBROUTINE  MAT__DGEMM_C_EQ_AXB_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      ROW,SUM,COL,
     +                      A,B,
     +
     +                              C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__DGEMM_C_EQ_AXB_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : Highly optimized matrix multiplication routine
C                performing the matrix multiplication:
C
C                                C = A * B
C
C                Matrices A and B need not be distinct, but matrix C
C                must be distinct from A and B and all matrices must
C                be of floating point type.
C
C                Code derived from the 'dgemm' Level 3 Blas routine.
C
C                To tune this routine for maximum performance, pass
C                ROWB,SUMB and COLB in argument and test for several
C                sizes of these variables.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    DDROWA,DDCOLA
         INTEGER    DDROWB,DDCOLB
         INTEGER    DDROWC,DDCOLC
         INTEGER    I,J,K
         INTEGER    I1,I2,J1,K1,II1,II2,JJ1,KK1
         INTEGER    IDEPTH,JDEPTH
         INTEGER    II,JJ,KK
         INTEGER    IIEND1,IIEND2,JJEND1,JJEND2,KKEND
         INTEGER    ILEN,JLEN
         INTEGER    ISPAN,JSPAN,KSPAN
         INTEGER    ROW,SUM,COL
         INTEGER    ROWB,SUMB,COLB

C Original value for the parameter (supposed to be for the IBM RS/6000)
C         PARAMETER  (ROWB=64,COLB=64,SUMB=64)
C Parameter for the 256kB cache on the IBM RS/6000 models 590 and 990
         PARAMETER  (ROWB=96,COLB=96,SUMB=96)
C Parameter for the 16kB cache on the Sparc
C         PARAMETER   (ROWB=24,COLB=24,SUMB=24)
C Parameter for the 4MB L2 cache on the SGI POWER CHALLANGE
C         PARAMETER  (ROWB=384,COLB=384,SUMB=384)
C Parameter for the 8kB cache on the DEC Alpha
C         PARAMETER  (ROWB=16,COLB=16,SUMB=16)
C Optimized parameters for the QTP crunch
C         PARAMETER   (ROWB=46,COLB=46,SUMB=46)

         DOUBLE PRECISION   E,F,G,H
         DOUBLE PRECISION   T11,T12,T21,T22
         DOUBLE PRECISION   ZERO

         DOUBLE PRECISION   A (1:DDROWA,DDCOLA)
         DOUBLE PRECISION   B (1:DDROWB,DDCOLB)
         DOUBLE PRECISION   C (1:DDROWC,DDCOLC)

         DOUBLE PRECISION   ASUB (1:SUMB,1:ROWB)

         PARAMETER   (ZERO = 0.D0)
C
C
C-----------------------------------------------------------------------
C
C
C             ...proceed.
C
C
         DO J = 1,COL
         DO I = 1,ROW
            C (I,J) = ZERO
         END DO
         END DO

         DO 2000 K = 1,SUM,SUMB
            K1 = K - 1
            KSPAN = MIN0 (SUMB,SUM-K1)
            KKEND = K1 + KSPAN

            DO 1000 I = 1,ROW,ROWB
               I1 = I - 1
               I2 = I - 2
               IDEPTH = 2
               ISPAN = MIN0 (ROWB,ROW-I1)
               IIEND2 = I1 + ISPAN
               ILEN = IDEPTH * (ISPAN/IDEPTH)
               IIEND1 = I1 + ILEN
               DO II = I,IIEND2
                  II1 = II - I1
                  DO KK = K,KKEND
                     ASUB (KK-K1,II1) = A (II,KK)
                  END DO
               END DO

               DO 3000 J = 1,COL,COLB
                  J1 = J - 1
                  JDEPTH = 2
                  JSPAN = MIN0 (COLB,COL-J1)
                  JJEND2 = J1 + JSPAN
                  JLEN = JDEPTH * (JSPAN/JDEPTH)
                  JJEND1 = J1 + JLEN

                  DO JJ = J,JJEND1,JDEPTH
                     JJ1 = JJ + 1
                     DO II = I,IIEND1,IDEPTH
                        II1 = II - I1
                        II2 = II - I2
                        T11 = ZERO
                        T21 = ZERO
                        T12 = ZERO
                        T22 = ZERO
                        DO KK = K,KKEND
                           KK1 = KK - K1
                           E = ASUB (KK1,II1)
                           F = ASUB (KK1,II2)
                           G = B (KK,JJ)
                           H = B (KK,JJ1)
                           T11 = T11 + E * G
                           T21 = T21 + F * G
                           T12 = T12 + E * H
                           T22 = T22 + F * H
                        END DO
                        C(II,JJ) = C(II,JJ) + T11
                        C(II+1,JJ) = C(II+1,JJ) + T21
                        C(II,JJ1) = C(II,JJ1) + T12
                        C(II+1,JJ1) = C(II+1,JJ1) + T22
                     END DO
                     IF (ILEN.LT.ISPAN) THEN
                         DO II = IIEND1+1,IIEND2
                            II1 = II - I1
                            T11 = ZERO
                            T12 = ZERO
                            DO KK = K,KKEND
                               E = ASUB (KK-K1,II1)
                               T11 = T11 + E * B (KK,JJ)
                               T12 = T12 + E * B (KK,JJ1)
                            END DO
                            C(II,JJ) = C(II,JJ) + T11
                            C(II,JJ1) = C(II,JJ1) + T12
                         END DO
                     END IF
                  END DO

                  IF (JLEN.LT.JSPAN) THEN
                      DO JJ = JJEND1+1,JJEND2
                         DO II = I,IIEND1,IDEPTH
                            II1 = II - I1
                            II2 = II - I2
                            T11 = ZERO
                            T21 = ZERO
                            DO KK = K,KKEND
                               KK1 = KK - K1
                               G = B (KK,JJ)
                               T11 = T11 + ASUB (KK1,II1) * G
                               T21 = T21 + ASUB (KK1,II2) * G
                            END DO
                            C(II,JJ) = C(II,JJ) + T11
                            C(II+1,JJ) = C(II+1,JJ) + T21
                         END DO
                         IF (ILEN.LT.ISPAN) THEN
                             DO II = IIEND1+1,IIEND2
                                II1 = II - I1
                                T11 = ZERO
                                DO KK = K,KKEND
                                   T11 = T11 + ASUB(KK-K1,II1)*B(KK,JJ)
                                END DO
                                C(II,JJ) = C(II,JJ) + T11
                             END DO
                         END IF
                      END DO
                  END IF

 3000          CONTINUE
 1000       CONTINUE
 2000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
