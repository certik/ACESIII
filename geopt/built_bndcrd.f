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
      SUBROUTINE BULT_BNDCRD(CARTCOORD, BMATRX, DISTAB, ICON1,  
     &                       ICON2, IBNDS, TOTREDNCO, NRATMS)
C         
C Setup the bond stretching B-matrix elements. 
C B(*,i,j)(x,y,z) = (x_i,y_i,z_i - x_j,y_j,z_j)/|Rij|
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER TOTREDNCO
      DIMENSION CARTCOORD(3*NRATMS), BMATRX(TOTREDNCO, 3*NRATMS)
C     
      DISTAB = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
C
      DO 10 IXYZ = 1, 3
C
         DIFFAB = CARTCOORD((3*ICON1 - 3) + IXYZ)
     &             - CARTCOORD((3*ICON2 -3) + IXYZ)
         BMATRX(IBNDS, (3*ICON1 - 3) + IXYZ) =  DIFFAB/DISTAB
         BMATRX(IBNDS, (3*ICON2 - 3) + IXYZ) = -DIFFAB/DISTAB
C
 10   CONTINUE

      RETURN
      END
