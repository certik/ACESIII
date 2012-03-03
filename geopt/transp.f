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

      SUBROUTINE TRANSP(A,B,NUM,DIS)
C
C    THIS ROUTINE PERFORMS A GENERAL MATRIX  TRANSPOSITION 
C
C       B(P,Q) <---  A(Q,P)
C
C    INPUT:
C             A  ......   MATRIX TO BE TRANSPOSED
C             NUM .....   LENGTH OF SECOND INDEX OF A
C                         (NUMBER OF DISTRIBUTIONS IN A)
C             DIS .....   LENGTH OF FIRST INDEX OF A
C                         (DISTRIBUTION LENGTH IN A)
C
C    OUTPUT:
C             B .......   TRANSPOSED MATRIX
C
CEND
C
C  CODED JULY/90  JG
C  unrolled 8/95 SG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS 
      DIMENSION A(DIS,NUM),B(NUM,DIS)
C
      IST = MOD(DIS,8)
      DO 1 IAB = 1,IST
        DO 2 IJ = 1,NUM
          B(IJ,IAB)=A(IAB,IJ)
 2      CONTINUE
 1    CONTINUE
      DO 10 IAB=IST+1,DIS,8
        DO 20 IJ=1,NUM
          B(IJ,IAB)=A(IAB,IJ)
          B(IJ,IAB+1)=A(IAB+1,IJ)
          B(IJ,IAB+2)=A(IAB+2,IJ)
          B(IJ,IAB+3)=A(IAB+3,IJ)
          B(IJ,IAB+4)=A(IAB+4,IJ)
          B(IJ,IAB+5)=A(IAB+5,IJ)
          B(IJ,IAB+6)=A(IAB+6,IJ)
          B(IJ,IAB+7)=A(IAB+7,IJ)
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
