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
      INTEGER FUNCTION IFINDNE(N,A,ISKIP,TARGET,TOLER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
      IOFF=1
      IPASS=1
c1     DIFF=ABS(A(IOFF)-a(max(1,ioff-1)))
1     DIFF=ABS(A(IOFF)-TARGET)
      IF(DIFF.GT.TOLER.OR.IPASS.EQ.N+1)THEN
       IFINDNE=IPASS
      ELSE
       IPASS=IPASS+1
       IOFF=IOFF+ISKIP
       GOTO 1
      ENDIF
      RETURN
      END
