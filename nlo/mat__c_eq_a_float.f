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
         SUBROUTINE  MAT__C_EQ_A_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWC,DDCOLC,
     +                      ROW,COL,
     +                      A,
     +
     +                              C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_A_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation copies a two-dimensional matrix A of
C                floating point type to matrix C.
C
C                The dimensions of both matrices need not be the same.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDROWA,DDCOLA
         INTEGER    DDROWC,DDCOLC
         INTEGER    I,J
         INTEGER    ROW,COL

         DOUBLE PRECISION   A (1:DDROWA,1:DDCOLA)
         DOUBLE PRECISION   C (1:DDROWC,1:DDCOLC)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (ROW.GT.DDROWA  .OR.  COL.GT.DDCOLA)  THEN
              WRITE (1,*) ' Dimensions of matrix A too small: '
              WRITE (1,*) ' mat__c_eq_a_float '
              WRITE (1,*) ' DDROWA,DDCOLA,ROW,COL = ',
     +                      DDROWA,DDCOLA,ROW,COL
              WRITE (*,*) ' Dimensions of matrix A too small: '
              WRITE (*,*) ' mat__c_eq_a_float '
              WRITE (*,*) ' DDROWA,DDCOLA,ROW,COL = ',
     +                      DDROWA,DDCOLA,ROW,COL
              STOP
         END IF

         IF  (ROW.GT.DDROWC  .OR.  COL.GT.DDCOLC)  THEN
              WRITE (1,*) ' Dimensions of matrix C too small: '
              WRITE (1,*) ' mat__c_eq_a_float '
              WRITE (1,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                      DDROWC,DDCOLC,ROW,COL
              WRITE (*,*) ' Dimensions of matrix C too small: '
              WRITE (*,*) ' mat__c_eq_a_float '
              WRITE (*,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                      DDROWC,DDCOLC,ROW,COL
              STOP
         END IF
C
C
C             ...loop over all matrix elements.
C
C
         DO J = 1,COL
         DO I = 1,ROW  
            C (I,J) = A (I,J)
         END DO
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
