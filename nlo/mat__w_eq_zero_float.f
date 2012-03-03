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
         SUBROUTINE  MAT__W_EQ_ZERO_FLOAT
     +
     +                    ( DDVEC,
     +                      VEC,
     +
     +                              W )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__W_EQ_ZERO_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation builds a zero vector W of floating
C                point type.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDVEC
         INTEGER    I
         INTEGER    VEC

         DOUBLE PRECISION   W (1:DDVEC)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimension.
C
C
         IF  ( VEC .GT. DDVEC )  THEN

               WRITE (1,*) ' Dimension of vector W too small: '
               WRITE (1,*) ' mat__w_eq_zero_float '
               WRITE (1,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               WRITE (*,*) ' Dimension of vector W too small: '
               WRITE (*,*) ' mat__w_eq_zero_float '
               WRITE (*,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               STOP

         END IF
C
C
C             ...build zero vector.
C
C
         DO  10  I = 1,VEC  
             W (I) = 0.D0
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
