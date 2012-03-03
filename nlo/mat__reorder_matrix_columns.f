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
         SUBROUTINE  MAT__REORDER_MATRIX_COLUMNS
     +
     +                    ( DDROW,DDCOL,
     +                      DDORDER,
     +                      DDUSED,
     +                      DDWORK,
     +                      ROW,COL,
     +                      ORDER,
     +                      USED,
     +                      WORK,
     +                      C )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__REORDER_MATRIX_COLUMNS
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This routine reorders columns of input matrix C
C                according to a new ordering defined in array ORDER.
C                The ordering is done in place, hence C is overwritten
C                by the reordered matrix. Definition of the ordering
C                array is in the active sense:
C
C                      ORDER (old column #) = new column #
C
C                Example:
C
C                Let C be:
C
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C                          1  2  3  4  5  6  7  8
C
C                then, if ORDER is equal to:
C
C                              ORDER (1) = 5
C                              ORDER (2) = 2
C                              ORDER (3) = 1
C                              ORDER (4) = 6
C                              ORDER (5) = 8
C                              ORDER (6) = 3
C                              ORDER (7) = 7
C                              ORDER (8) = 4
C
C                the matrix C on exit looks like:
C
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
C                          3  2  6  8  1  4  7  5
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

         INTEGER     DDORDER
         INTEGER     DDROW,DDCOL
         INTEGER     DDUSED
         INTEGER     DDWORK
         INTEGER     I,J,N
         INTEGER     NEWCOL
         INTEGER     ROW,COL

         INTEGER     ORDER (1:DDORDER)
         INTEGER     USED  (1:DDUSED)

         DOUBLE PRECISION  CVAL

         DOUBLE PRECISION  WORK (1:DDWORK)

         DOUBLE PRECISION  C (1:DDROW,1:DDCOL)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (ROW.GT.DDROW .OR. COL.GT.DDCOL) THEN
              WRITE (*,*) ' Dimensions of matrix C too small! '
              WRITE (*,*) ' mat__reorder_matrix_columns '
              WRITE (*,*) ' DDROW,DDCOL,ROW,COL = ',DDROW,DDCOL,ROW,COL
              WRITE (1,*) ' Dimensions of matrix C too small! '
              WRITE (1,*) ' mat__reorder_matrix_columns '
              WRITE (1,*) ' DDROW,DDCOL,ROW,COL = ',DDROW,DDCOL,ROW,COL
              STOP
         END IF

         IF  (ROW.GT.DDWORK) THEN
              WRITE (*,*) ' Dimensions of array WORK too small! '
              WRITE (*,*) ' mat__reorder_matrix_columns '
              WRITE (*,*) ' DDWORK,ROW = ',DDWORK,ROW
              WRITE (1,*) ' Dimensions of array WORK too small! '
              WRITE (1,*) ' mat__reorder_matrix_columns '
              WRITE (1,*) ' DDWORK,ROW = ',DDWORK,ROW
              STOP
         END IF

         IF  (COL.GT.DDORDER) THEN
              WRITE (*,*) ' Dimensions of array ORDER too small! '
              WRITE (*,*) ' mat__reorder_matrix_columns '
              WRITE (*,*) ' DDORDER,COL = ',DDORDER,COL
              WRITE (1,*) ' Dimensions of array ORDER too small! '
              WRITE (1,*) ' mat__reorder_matrix_columns '
              WRITE (1,*) ' DDORDER,COL = ',DDORDER,COL
              STOP
         END IF

         IF  (COL.GT.DDUSED) THEN
              WRITE (*,*) ' Dimensions of array USED too small! '
              WRITE (*,*) ' mat__reorder_matrix_columns '
              WRITE (*,*) ' DDUSED,COL = ',DDUSED,COL
              WRITE (1,*) ' Dimensions of array USED too small! '
              WRITE (1,*) ' mat__reorder_matrix_columns '
              WRITE (1,*) ' DDUSED,COL = ',DDUSED,COL
              STOP
         END IF
C
C
C             ...check the reordering array passed. It might have
C                column indices out of range, that is numbers > COL.
C
C
         DO J = 1,COL
            NEWCOL = ORDER (J)
            IF (NEWCOL.GT.COL) THEN
                WRITE (*,*) ' Column # out of meaningful range! '
                WRITE (*,*) ' mat__reorder_matrix_columns '
                WRITE (*,*) ' NEWCOL,J,COL = ',NEWCOL,J,COL
                WRITE (1,*) ' Column # out of meaningful range! '
                WRITE (1,*) ' mat__reorder_matrix_columns '
                WRITE (1,*) ' NEWCOL,J,COL = ',NEWCOL,J,COL
                STOP
            END IF
         END DO
C
C
C             ...initialize used column array and start loop over
C                all cycles.
C
C
         DO J = 1,COL
            USED (J) = 0
         END DO

         DO 1000 J = 1,COL
            IF (USED (J).EQ.0) THEN
                NEWCOL = ORDER (J)
C
C
C             ...handle current cycle, if not of length 1.
C                Inner loop over dummy index N runs over all
C                # of columns to ensure that end of cycle is
C                reached.
C
C
                IF (NEWCOL.NE.J) THEN

                    USED (NEWCOL) = 1
                    DO I = 1,ROW
                       WORK (I) = C (I,NEWCOL)
                       C (I,NEWCOL) = C (I,J)
                    END DO

                    DO N = 1,COL
                       NEWCOL = ORDER (NEWCOL)
                       USED (NEWCOL) = 1
                       DO I = 1,ROW
                          CVAL = C (I,NEWCOL)
                          C (I,NEWCOL) = WORK (I)
                          WORK (I) = CVAL
                       END DO
                       IF (NEWCOL.EQ.J) GOTO 1000
                    END DO

                END IF

            END IF
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
