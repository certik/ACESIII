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
      SUBROUTINE INITP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,EOR,DSTRT
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
      PARAMETER (LENK  = 8*MXCBF + 9 + 7 + 8 + 3 + 7*28 + 64 + 8*8*2
     &                 + 8*2 + 100 + MXTNSH*8) 
C-----------------------------------------------------------------------
      COMMON /INDX/ PC(512),K(LENK)
c     COMMON /INDX/ PC(512),K(5011)
c     COMMON /INDX/ PC(512),K(512)
C     INITIALIZE THE POPULATION COUNT FACTOR
      K(1)=0
      N=1
      DO 1 I=1,9
        DO 2 J=1,N
    2   K(J+N) = K(J)+1
    1 N = N + N
      DO 3 I=1,512
    3 PC(I) = DFLOAT(1 - 2*AND(K(I),1))
      RETURN
      END
