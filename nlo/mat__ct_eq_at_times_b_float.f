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
         SUBROUTINE  MAT__CT_EQ_AT_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      ROW,SUM,COL,
     +                      A,B,
     +
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__CT_EQ_AT_TIMES_B_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation multiplies the transpose of the
C                two-dimensional matrix A by the two-dimensional matrix
C                B. The result will be placed in transposed form into C.
C                Matrices A and B need not be distinct, but matrix C
C                must be distinct from A and B and all matrices must
C                be of floating point type.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    DDROWA,DDCOLA
         INTEGER    DDROWB,DDCOLB
         INTEGER    DDROWC,DDCOLC
         INTEGER    ROW,SUM,COL

         DOUBLE PRECISION   A (1:DDROWA,1:DDCOLA)
         DOUBLE PRECISION   B (1:DDROWB,1:DDCOLB)
         DOUBLE PRECISION   C (1:DDROWC,1:DDCOLC)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (ROW.GT.DDCOLA .OR. SUM.GT.DDROWA)  THEN
              WRITE (1,*) ' Dimensions of matrix A too small: '
              WRITE (1,*) ' mat__ct_eq_at_times_b_float '
              WRITE (1,*) ' DDCOLA,DDROWA,ROW,SUM = ',
     +                      DDCOLA,DDROWA,ROW,SUM
              WRITE (*,*) ' Dimensions of matrix A too small: '
              WRITE (*,*) ' mat__ct_eq_at_times_b_float '
              WRITE (*,*) ' DDCOLA,DDROWA,ROW,SUM = ',
     +                      DDCOLA,DDROWA,ROW,SUM
              STOP
         END IF

         IF  (SUM.GT.DDROWB .OR. COL.GT.DDCOLB)  THEN
              WRITE (1,*) ' Dimensions of matrix B too small: '
              WRITE (1,*) ' mat__ct_eq_at_times_b_float '
              WRITE (1,*) ' DDROWB,DDCOLB,SUM,COL = ',
     +                      DDROWB,DDCOLB,SUM,COL
              WRITE (*,*) ' Dimensions of matrix B too small: '
              WRITE (*,*) ' mat__ct_eq_at_times_b_float '
              WRITE (*,*) ' DDROWB,DDCOLB,SUM,COL = ',
     +                      DDROWB,DDCOLB,SUM,COL
              STOP
         END IF

         IF  (ROW.GT.DDCOLC .OR. COL.GT.DDROWC)  THEN
              WRITE (1,*) ' Dimensions of matrix C too small: '
              WRITE (1,*) ' mat__ct_eq_at_times_b_float '
              WRITE (1,*) ' DDCOLC,DDROWC,ROW,COL = ',
     +                      DDCOLC,DDROWC,ROW,COL
              WRITE (*,*) ' Dimensions of matrix C too small: '
              WRITE (*,*) ' mat__ct_eq_at_times_b_float '
              WRITE (*,*) ' DDCOLC,DDROWC,ROW,COL = ',
     +                      DDCOLC,DDROWC,ROW,COL
              STOP
         END IF
C
C
C             ...highly optimized routine.
C
C
         CALL  MAT__DGEMM_CT_EQ_ATXB_FLOAT
     +
     +              ( DDROWA,DDCOLA,
     +                DDROWB,DDCOLB,
     +                DDROWC,DDCOLC,
     +                ROW,SUM,COL,
     +                A,B,
     +
     +                       C )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
