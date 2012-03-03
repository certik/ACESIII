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
      SUBROUTINE SWPORD(NLIST1,NLIST2,VECIN,VECOUT,NSIZ)
C 
C SUBROUTINE SWAPS ORDERING OF COORDINATES IN VECIN DEPENDING ON
C  POINTER LISTS NLIST1 AND NLIST2.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NLIST1(NSIZ),NLIST2(NSIZ),VECIN(3*NSIZ),VECOUT(3*NSIZ)
      CALL ZERO(VECOUT,3*NSIZ)
      DO 10 I=1,NSIZ
       IP1=3*NLIST1(I)-2
       IP2=3*NLIST2(I)-2
       CALL VADD(VECOUT(IP1),VECOUT(IP1),VECIN(IP2),3,1.D0)
 10   CONTINUE
      RETURN
      END
