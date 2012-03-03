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
C
      SUBROUTINE CSHIFTPRI(CSHIFT, INDEX)
C
C Print out the total shielding tensor and orientation of the shielding
C tensor. Coded by Ajith 04/96.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION CSHIFT (3, 3)
C
      CHARACTER*5 COORD(3)
C
#include "mxatms.par"
#include "cbchar.com"
#include "coord.com"
C
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
C
      DATA COORD/'x    ','y    ','z    '/
C
      WRITE (6, 10)
C
      DO 5 I = 1, 3
         WRITE(6, 20) ZSYM(INDEX), COORD(I), (CSHIFT(I,J),J=1,3)
    5 CONTINUE
C
 10   FORMAT(/)
 20   FORMAT(15X, A2, 1X, A2, 3F12.6)
C
      RETURN
      END
