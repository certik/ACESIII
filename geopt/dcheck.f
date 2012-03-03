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
      SUBROUTINE DCHECK(RREF,RVEC,ISIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RVEC(ISIZE),RREF(ISIZE)
      PI=DACOS(-1.D0)
      DO 55 I=6,ISIZE,3
       IF(DABS(RVEC(I)).GT.3.D0)THEN
        IF(SIGN(1.D0,RVEC(I)).NE.SIGN(1.D0,RREF(I)))
     &  RVEC(I)=SIGN(1.D0,RREF(I))*(2.D0*PI-DABS(RVEC(I)))
       ENDIF
 55   CONTINUE
      RETURN
      END
