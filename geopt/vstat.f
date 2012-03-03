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
      SUBROUTINE VSTAT(V,ZQ,LENGTH)
C
C RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
C     ZQ(1)  Largest absolute magnitude
C     ZQ(2)  Smallest absolute magnitude
C     ZQ(3)  Largest value
C     ZQ(4)  Smallest value
C     ZQ(5)  2-norm
C     ZQ(6)  Dynamic range of the vector (abs. min. - abs. max.)
C   

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LENGTH),ZQ(6)
      U=0.D0
      CALL ZERO(ZQ,6)
      ZQ(2)=DABS(V(1))
      ZQ(4)=V(1)
      DO 20 I=1,LENGTH
      ZQ(1)=MAX(ZQ(1),DABS(V(I)))
      ZQ(2)=MIN(ZQ(2),DABS(V(I)))
      ZQ(3)=MAX(ZQ(3),V(I))
      ZQ(4)=MIN(ZQ(4),V(I))
20    U=U+V(I)*V(I)
      If (Length .ne. 0) ZQ(5)=DSQRT(U/LENGTH)
      ZQ(6)=ZQ(2)-ZQ(1)
      RETURN
      END
