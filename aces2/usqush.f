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
      SUBROUTINE USQUSH(V,ISIZE)
C
C "UNSQUASHES" A VECTOR OF LENGTH ISIZE.
C
      DOUBLE PRECISION V(ISIZE+6)
      DO 10 I=ISIZE,4,-1
      V(I+6)=V(I)
 10   CONTINUE
      V(8)=V(3)
      V(7)=V(2)
      V(4)=V(1)
      V(9)=0.D0
      V(6)=0.D0
      V(5)=0.D0
      DO 20 I=1,3
      V(I)=0.D0
 20   CONTINUE
      RETURN
      END
