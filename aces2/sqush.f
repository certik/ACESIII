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
      SUBROUTINE SQUSH(VEC,LEN)
C
C SUBROUTINE USED TO "SQUASH" R VECTOR.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(LEN)
      VEC(1)=VEC(4)
      VEC(2)=VEC(7)
      VEC(3)=VEC(8)
      DO 10 I=4,LEN-6
      VEC(I)=VEC(I+6)
 10   CONTINUE
      CALL ZERO(VEC(MAX(LEN-5,2)),6)
      RETURN
      END
