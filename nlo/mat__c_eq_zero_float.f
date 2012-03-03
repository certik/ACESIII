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
         SUBROUTINE  MAT__C_EQ_ZERO_FLOAT
     +
     +                    ( DDROW,DDCOL,
     +                      ROW,COL,
     +
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_ZERO_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation builds a two-dimensional zero matrix C
C                of floating point type.
C
C                The zero matrix need not to be a square matrix.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    DDROW,DDCOL
         INTEGER    I,J
         INTEGER    ROW,COL

         DOUBLE PRECISION   ZERO

         DOUBLE PRECISION   C (1:DDROW,1:DDCOL)

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (ROW.GT.DDROW  .OR.  COL.GT.DDCOL)  THEN
              WRITE (1,*) ' Dimensions of matrix C too small: '
              WRITE (1,*) ' mat__c_eq_zero_float '
              WRITE (1,*) ' DDROW,DDCOL,ROW,COL = ',
     +                      DDROW,DDCOL,ROW,COL
              WRITE (*,*) ' Dimensions of matrix C too small: '
              WRITE (*,*) ' mat__c_eq_zero_float '
              WRITE (*,*) ' DDROW,DDCOL,ROW,COL = ',
     +                      DDROW,DDCOL,ROW,COL
              STOP
         END IF
C
C
C             ...build zero matrix.
C
C
         DO J = 1,COL
         DO I = 1,ROW  
            C (I,J) = ZERO
         END DO
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
