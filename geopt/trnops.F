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

C TRANSFORMS SYMMETRY OPERATIONS TO AN ORIENTATION SPECIFIED BY
C MATRIX "GOOFY".

      SUBROUTINE TRNOPS(OPS,GOOFY,IORDGP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION OPS(9*IORDGP),GOOFY(9),SCR(9)

#ifdef _ASSERT
      if (iordgp.lt.0) then
         print *, '@TRNOPS: Assertion failed.'
         print *, '         iordgp = ',iordgp
         call errex
      end if
#endif /* _ASSERT */

      DO I = 1, IORDGP
         IADD = -8 + (9*I)
         CALL XGEMM('N','N',3,3,3,
     &              1.0d0,GOOFY,    3,
     &                    OPS(IADD),3,
     &              0.0d0,SCR,      3)
         CALL XGEMM('N','T',3,3,3,
     &              1.0d0,SCR,      3,
     &                    GOOFY,    3,
     &              0.0d0,OPS(IADD),3)
      END DO

      RETURN
      END 

