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
*DECK ELTRAN
      SUBROUTINE ELTRAN (NM, N, LOW, IGH, A, INT, Z)
C***BEGIN PROLOGUE  ELTRAN
C***PURPOSE  Accumulates the stabilized elementary similarity
C            transformations used in the reduction of a real general
C            matrix to upper Hessenberg form by ELMHES.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      DOUBLE PRECISION (ELTRAN-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine accumulates the stabilized elementary
C     similarity transformations used in the reduction of a
C     REAL*8 GENERAL matrix to upper Hessenberg form by  ELMHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  BALANC.  If  BALANC  has not been
C          used, set LOW=1 and IGH equal to the order of the matrix, N.
C
C        A contains the multipliers which were used in the reduction
C          by  ELMHES  in its lower triangle below the subdiagonal.
C          A is a two-dimensional REAL*8 array, dimensioned A(NM,IGH).
C
C        INT contains information on the rows and columns interchanged
C          in the reduction by  ELMHES.  Only elements LOW through IGH
C          are used.  INT is a one-dimensional INTEGER array,
C          dimensioned INT(IGH).
C
C     On OUTPUT
C
C        Z contains the transformation matrix produced in the reduction
C          by  ELMHES.  Z is a two-dimensional REAL*8 array, dimensioned
C          Z(NM,N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  ELTRAN
C
#ifdef ACES64
      INTEGER*8 I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      INTEGER*8 INT(*)
#else
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      INTEGER INT(*)
#endif
      DOUBLE PRECISION A(NM,*),Z(NM,*)
C
C***FIRST EXECUTABLE STATEMENT  ELTRAN
      DO 80 J = 1, N
C
         DO 60 I = 1, N
   60    Z(I,J) = 0.0D0
C
         Z(J,J) = 1.0D0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0D0
  130    CONTINUE
C
         Z(I,MP) = 1.0D0
  140 CONTINUE
C
  200 RETURN
      END
