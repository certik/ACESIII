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

C THIS ROUTINE NEGATES THE FIRST LEN ELEMENTS OF A VECTOR V.

      SUBROUTINE VMINUS(V,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LEN)
      if (len.lt.1) return
      DO I = 1, LEN
         V(I) = -V(I)
      END DO
      RETURN
      END
