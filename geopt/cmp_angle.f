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
      SUBROUTINE CMP_ANGLE(ILFTEND, IMIDLE, IRHTEND, CARTCOORD,
     &                     VLUANGLE, NRATMS) 
C
C This is a simple routine to calculate the angle between two vectors,
C ILFTEND and IRHTEND and the center is IMIDLE. The angle is retrun
C in VLUANGLE.
C      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION CARTCOORD(3*NRATMS), VEC1(3), VEC2(3)


      CALL VEC(CARTCOORD(IMIDLE*3 - 2), CARTCOORD(ILFTEND*3 - 2), 
     &         VEC1, 1)
C
      CALL VEC(CARTCOORD(IMIDLE*3 - 2), CARTCOORD(IRHTEND*3 - 2), 
     &         VEC2, 1)

      VLUANGLE = ANGLE(VEC1, VEC2, 3)
C
      RETURN
      END
