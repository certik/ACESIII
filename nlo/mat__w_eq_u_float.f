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
         SUBROUTINE  MAT__W_EQ_U_FLOAT
     +
     +                    ( DDVECU,
     +                      DDVECW,
     +                      VEC,
     +                      U,
     +
     +                              W )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__W_EQ_U_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation copies a vector U of floating point
C                type to vector W.
C
C                The dimensions of both vectors need not be the same.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDVECU
         INTEGER    DDVECW
         INTEGER    I
         INTEGER    VEC

         DOUBLE PRECISION   U (1:DDVECU)
         DOUBLE PRECISION   W (1:DDVECW)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  ( VEC .GT. DDVECU )  THEN

               WRITE (1,*) ' Dimension of vector U too small: '
               WRITE (1,*) ' mat__w_eq_u_float '
               WRITE (1,*) ' DDVECU,VEC = ',
     +                       DDVECU,VEC

               WRITE (*,*) ' Dimension of vector U too small: '
               WRITE (*,*) ' mat__w_eq_u_float '
               WRITE (*,*) ' DDVECU,VEC = ',
     +                       DDVECU,VEC

               STOP

         END IF

         IF  ( VEC .GT. DDVECW )  THEN

               WRITE (1,*) ' Dimension of vector W too small: '
               WRITE (1,*) ' mat__w_eq_u_float '
               WRITE (1,*) ' DDVECW,VEC = ',
     +                       DDVECW,VEC

               WRITE (*,*) ' Dimension of vector W too small: '
               WRITE (*,*) ' mat__w_eq_u_float '
               WRITE (*,*) ' DDVECW,VEC = ',
     +                       DDVECW,VEC

               STOP

         END IF
C
C
C             ...loop over all vector elements.
C
C
         DO  10  I = 1,VEC  
             W (I) = U (I)
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
