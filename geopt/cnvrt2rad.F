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
      SUBROUTINE CNVRT2RAD(R,NX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      DIMENSION R(NX)
      ATOR=DACOS(-1.D0)/180.D0
      ATOB=0.529177249D0
      DO IX = 8, NX-1, 3
         IF (IX.NE.8) R(IX+1) = R(IX+1)*ATOR
         R(IX) = R(IX)*ATOR
      END DO
      IF (IFLAGS(78).EQ.0) THEN
         DO IX = 4, NX-2, 3
            R(IX) = R(IX)/ATOB
         END DO
      END IF
      RETURN
      END
