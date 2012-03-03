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

C This routine expands a triangularly packed vector of numbers
C into a square matrix.
C
C  WPACK((NDIM*(NDIM+1))/2) ==> WFULL(NDIM,NDIM)

      SUBROUTINE EXPND2(WPACK,WFULL,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WPACK((NDIM*(NDIM+1))/2),WFULL(NDIM,NDIM)
#ifdef _ASSERT
      if (ndim.lt.0) then
         print *, '@EXPND2: Assertion failed.'
         print *, '         ndim = ',ndim
         call errex
      end if
#endif /* _ASSERT */
      ITHRU = 0
      DO I = 1, NDIM
         DO J = 1, I
            ITHRU = ITHRU + 1
            WFULL(I,J) = WPACK(ITHRU)
            WFULL(J,I) = WPACK(ITHRU)
         END DO
      END DO
      RETURN
      END
