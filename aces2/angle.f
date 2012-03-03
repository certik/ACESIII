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
      DOUBLE PRECISION FUNCTION ANGLE(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N)
      DTOR=DACOS(-1.D0)/180.D0
      Z = ddot(N,A,1,B,1)
C
C MAKE SURE A VALUE OF SLIGHTLY MORE THAN ONE DOES NOT GET PASSED
C  TO DACOS
C
      ANGLE=DACOSX(Z,1.D-10)/DTOR
      RETURN
      END
