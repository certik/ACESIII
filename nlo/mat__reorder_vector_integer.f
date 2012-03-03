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
         SUBROUTINE  MAT__REORDER_VECTOR_INTEGER
     +
     +                    ( DDVEC,
     +                      DDORDER,
     +                      DDWORK,
     +                      VEC,
     +                      ORDER,
     +                      WORK,
     +                      W )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__REORDER_VECTOR_INTEGER
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This routine reorders the elements of input vector W
C                of integer type according to a new ordering defined
C                in array ORDER. Definition of the ordering array is in
C                the active sense:
C
C                        ORDER (old element) = new element
C
C                Example:
C
C                Let W be:
C
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
C                the vector W on exit looks like:
C
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
         INTEGER     DDVEC
         INTEGER     DDWORK
         INTEGER     I
         INTEGER     VEC

         INTEGER     ORDER (1:DDORDER)
         INTEGER     W     (1:DDVEC)
         INTEGER     WORK  (1:DDWORK)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  (VEC.GT.DDVEC) THEN
              WRITE (*,*) ' Dimensions of vector W too small! '
              WRITE (*,*) ' nlo__reorder_vector_integer '
              WRITE (*,*) ' DDVEC,VEC = ',DDVEC,VEC
              WRITE (1,*) ' Dimensions of vector W too small! '
              WRITE (1,*) ' nlo__reorder_vector_integer '
              WRITE (1,*) ' DDVEC,VEC = ',DDVEC,VEC
              STOP
         END IF

         IF  (VEC.GT.DDWORK) THEN
              WRITE (*,*) ' Dimensions of array WORK too small! '
              WRITE (*,*) ' nlo__reorder_vector_integer '
              WRITE (*,*) ' DDWORK,VEC = ',DDWORK,VEC
              WRITE (1,*) ' Dimensions of array WORK too small! '
              WRITE (1,*) ' nlo__reorder_vector_integer '
              WRITE (1,*) ' DDWORK,VEC = ',DDWORK,VEC
              STOP
         END IF

         IF  (VEC.GT.DDORDER) THEN
              WRITE (*,*) ' Dimensions of array ORDER too small! '
              WRITE (*,*) ' nlo__reorder_vector_integer '
              WRITE (*,*) ' DDORDER,VEC = ',DDORDER,VEC
              WRITE (1,*) ' Dimensions of array ORDER too small! '
              WRITE (1,*) ' nlo__reorder_vector_integer '
              WRITE (1,*) ' DDORDER,VEC = ',DDORDER,VEC
              STOP
         END IF
C
C
C             ...check the reordering array passed. It might have
C                element indices out of range, that is numbers > VEC.
C
C
         DO I = 1,VEC
            IF (ORDER (I).GT.VEC) THEN
                WRITE (*,*) ' Element index out of meaningful range! '
                WRITE (*,*) ' nlo__reorder_vector_integer '
                WRITE (*,*) ' I,Element idx I,VEC = ',I,ORDER (I),VEC
                WRITE (1,*) ' Element index out of meaningful range! '
                WRITE (1,*) ' nlo__reorder_vector_integer '
                WRITE (1,*) ' I,Element idx I,VEC = ',I,ORDER (I),VEC
                STOP
            END IF
         END DO
C
C
C             ...reorder the vector elements in straightforward way.
C
C
         DO I = 1,VEC
            WORK (I) = W (I)
         END DO

         DO I = 1,VEC
            W (ORDER (I)) = WORK (I)
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
