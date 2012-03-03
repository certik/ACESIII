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
         SUBROUTINE  MAT__C_EQ_C_TRANSPOSED_FLOAT
     +
     +                    ( DDROWC,DDCOLC,
     +                      ROW,COL,
     +                      C )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_C_TRANSPOSED_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation in place transposes a 2-dimensional
C                matrix C of floating point type.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    DDROWC,DDCOLC
         INTEGER    I,J,K
         INTEGER    ROW,COL

         DOUBLE PRECISION   TEMP

         DOUBLE PRECISION   C (1:DDROWC,1:DDCOLC)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (ROW.GT.DDCOLC  .OR.  COL.GT.DDROWC)  THEN
              WRITE (1,*) ' Dimensions of matrix C too small: '
              WRITE (1,*) ' mat__c_eq_c_transposed_float '
              WRITE (1,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                      DDROWC,DDCOLC,ROW,COL
              WRITE (*,*) ' Dimensions of matrix C too small: '
              WRITE (*,*) ' mat__c_eq_c_transposed_float '
              WRITE (*,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                      DDROWC,DDCOLC,ROW,COL
              STOP
         END IF
C
C
C             ...transpose common upper left square.
C
C
         K = MIN0 (ROW,COL)

         DO I = 1,K
         DO J = I,K
            TEMP = C (I,J)
            C (I,J) = C (J,I)
            C (J,I) = TEMP
         END DO
         END DO
C
C
C             ...transpose rest of C.
C
C
         IF (ROW.GT.COL) THEN

             DO J = 1,COL
             DO I = COL+1,ROW
                C (J,I) = C (I,J)
             END DO
             END DO

         ELSE IF (ROW.LT.COL) THEN

             DO I = 1,ROW
             DO J = ROW+1,COL
                C (J,I) = C (I,J)
             END DO
             END DO

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
