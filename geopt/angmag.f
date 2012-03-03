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
        DOUBLE PRECISION FUNCTION ANGMAG(ORBIT,D,IORDER,IRETURN)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IRETURN=0
        ANGMAG=0.D0
        DTOR=DACOS(-1.D0)/180.D0
        ORDER=DFLOAT(IORDER)
        A1=180.D0/ORDER
        IF(DABS(A1-90.D0).LT.1.D-3)THEN
          ANGMAG=0.D0
        RETURN
        ENDIF
        X=ORBIT**2-0.25D0*D**2/DSIN(DTOR*A1)**2
        TOP=D/DTAN(DTOR*A1)
        IF(X.LT.-1.D-12)THEN
C
C THE ROTATION IS IMPOSSIBLE
C
         IRETURN=1
         RETURN
        ENDIF
        IF(DABS(X).LT.1.D-6)THEN
          ANGMAG=90.D0
        ELSE
        BOT=2.D0*DSQRT(X)
          ANGMAG=DATAN(TOP/BOT)/DTOR
        ENDIF
        RETURN
        END         
