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
         SUBROUTINE  MAT__C_EQ_C_TIMES_A_FLOAT
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWA,DDCOLA,
     +                      DDWORK1,DDWORK2,
     +                      ROW,SUM,COL,
     +                      WORK1,WORK2,
     +                      A,
     +
     +                              C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_C_TIMES_A_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation multiplies the two-dimensional matrix C
C                by the two-dimensional matrix A and stores the result
C                in matrix C. Matrix is C is thus overwritten:
C
C                                C = C * A
C
C                All matrices must be of floating point type.
C
C                Elementwise we can write the multiplication as:
C
C                          C(ij)  =  sum  C(ik) * A(kj)
C                                     k
C
C                where the index ranges are:
C
C                           i index    ->   1 to ROW
C                           j index    ->   1 to COL
C                           k index    ->   1 to SUM
C
C                Thus overwriting matrix C makes only sense, if its
C                column dimension declared is not smaller than the
C                # of columns of matrix A, i.e. we must have:
C
C                               DDCOLC >= COL
C
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
         INTEGER    DDROWC,DDCOLC
         INTEGER    DDWORK1,DDWORK2
         INTEGER    I,J,K
         INTEGER    ROW,SUM,COL

         DOUBLE PRECISION   X,ZERO

         DOUBLE PRECISION   WORK1 (1:DDWORK1)
         DOUBLE PRECISION   WORK2 (1:DDWORK2)

         DOUBLE PRECISION   A (1:DDROWA,1:DDCOLA)
         DOUBLE PRECISION   C (1:DDROWC,1:DDCOLC)

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions for matrix A.
C
C
         IF  (SUM.GT.DDROWA .OR. COL.GT.DDCOLA)  THEN
              WRITE (*,*) ' Dimensions of matrix A too small! '
              WRITE (*,*) ' mat__c_eq_c_times_a_float '
              WRITE (*,*) ' DDROWA,DDCOLA,SUM,COL = ',
     +                      DDROWA,DDCOLA,SUM,COL
              WRITE (1,*) ' Dimensions of matrix A too small! '
              WRITE (1,*) ' mat__c_eq_c_times_a_float '
              WRITE (1,*) ' DDROWA,DDCOLA,SUM,COL = ',
     +                      DDROWA,DDCOLA,SUM,COL
              STOP
         END IF
C
C
C             ...check dimensions for matrix C.
C
C
         IF  (ROW.GT.DDROWC .OR. MAX0 (SUM,COL).GT.DDCOLC)  THEN
              WRITE (*,*) ' Dimensions of matrix C too small! '
              WRITE (*,*) ' mat__c_eq_c_times_a_float '
              WRITE (*,*) ' DDROWC,DDCOLC,ROW,SUM,COL = ',
     +                      DDROWC,DDCOLC,ROW,SUM,COL
              WRITE (1,*) ' Dimensions of matrix C too small! '
              WRITE (1,*) ' mat__c_eq_c_times_a_float '
              WRITE (1,*) ' DDROWC,DDCOLC,ROW,SUM,COL = ',
     +                      DDROWC,DDCOLC,ROW,SUM,COL
              STOP
         END IF
C
C
C             ...check dimensions for working vectors.
C
C
         IF  (SUM.GT.DDWORK1 .OR. COL.GT.DDWORK2)  THEN
              WRITE (*,*) ' Dimensions of working arrays too small! '
              WRITE (*,*) ' mat__c_eq_c_times_a_float '
              WRITE (*,*) ' DDWORK1,DDWORK2,SUM,COL = ',
     +                      DDWORK1,DDWORK2,SUM,COL
              WRITE (1,*) ' Dimensions of working arrays too small! '
              WRITE (1,*) ' mat__c_eq_c_times_a_float '
              WRITE (1,*) ' DDWORK1,DDWORK2,SUM,COL = ',
     +                      DDWORK1,DDWORK2,SUM,COL
              STOP
         END IF
C
C
C             ...proceed in most cache effective way.
C
C
         DO I = 1,ROW

            DO K = 1,SUM
               WORK1 (K) = C (I,K)
            END DO

            DO J = 1,COL
               X = ZERO
               DO K = 1,SUM
                  X = X + WORK1 (K) * A (K,J)
               END DO
               WORK2 (J) = X
            END DO

            DO J = 1,COL
               C (I,J) = WORK2 (J)
            END DO

         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
