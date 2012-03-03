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
      SUBROUTINE STPTR(N,NORD1,NORD2,IPTR)
C
C ROUTINE TO CONSTRUCT POINTER LIST FOR FINITE DIFFERENCE CALCULATIONS.
C  IPTR(I)=J MEANS THAT THE SYMMETRY OPERATION MAPS ATOM I INTO ATOM J.
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION NORD1(N),NORD2(N),IPTR(N)
      DO 10 I=1,N
       J=NORD1(I)
       K=NORD2(I)
       IPTR(J)=K
10    CONTINUE
      RETURN
      END
