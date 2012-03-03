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
         SUBROUTINE  MAT__PRINT_V_INTEGER_NOZEROS
     +
     +                    ( UNITID,
     +                      TITLE,
     +                      DDVEC,
     +                      VEC,
     +                      V )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__PRINT_V_INTEGER_NOZEROS
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation prints a vector V to the unit specified
C                by UNITID in integer format I8 . Only values below
C                100000000 can be printed.
C
C                The vector is printed in one column.                    
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
         INTEGER    UNITID
         INTEGER    VEC
         INTEGER    ZERO

         INTEGER   V (1:DDVEC)

         CHARACTER*8     S
         CHARACTER*(*)   TITLE

         DATA  ZERO  /0/
C
C
C------------------------------------------------------------------------
C
C
C             ...printout title.
C
C
         WRITE (UNITID,6000)
         WRITE (UNITID, *  ) TITLE
         WRITE (UNITID,7000)
C
C
C             ...immediate return if length is zero.
C
C
         IF  ( VEC.EQ.0 )  RETURN
C
C
C             ...check dimension.
C
C
         IF  ( VEC .GT. DDVEC )  THEN

               WRITE (1,*) ' Dimension of vector V too small: '
               WRITE (1,*) ' mat__print_v_integer_nozeros '
               WRITE (1,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               WRITE (*,*) ' Dimensions of vector V too small: '
               WRITE (*,*) ' mat__print_v_integer_nozeros '
               WRITE (*,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               STOP

         END IF
C
C
C             ...print vector.
C
C
         DO  10  I = 1,VEC

             IF  ( V(I) .EQ. ZERO )  THEN
                   S = ' '
             ELSE 
                   WRITE (S,8000) V(I)
             END IF

             WRITE (UNITID,9000) I,S

   10    CONTINUE
C
C
C             ...formats of printing.
C
C
 6000    FORMAT  (/)
 7000    FORMAT  (/)
 8000    FORMAT  (I8)
 9000    FORMAT  (1X,I8,A8)
C
C
C             ...ready!
C
C
         RETURN
         END
