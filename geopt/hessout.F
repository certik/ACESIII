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
      SUBROUTINE HESSOUT(X,IR,IC,NSYM)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "coord.com"
C
C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
#include "cbchar.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
#include "io_units.par"
      DIMENSION X(IR,IC)
C
      JX=1
 1    ITOP=MIN(JX+4,IR)
      IF(NSYM.EQ.1)
     $   WRITE(LuOut,90)(VARNAM(ISQUASH(NOPTI(IX))),IX=JX,ITOP)
      IF(NSYM.EQ.0)WRITE(LuOut,90)(VARNAM(ISQUASH(IX)),IX=JX,ITOP)
 90   FORMAT(T15,A5,T25,A5,T35,A5,T45,A5,T55,A5)
      CALL OUTPUT(X, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#ifdef _NOSKIP
      DO 10 I=1,IR
         IF (NSYM .EQ. 1) WRITE (LuOut,44) VARNAM(ISQUASH(NOPTI(I))),
     $      (X(I,J),J=JX,MIN(IC,ITOP))
         IF (NSYM .EQ. 0) WRITE (LuOut,44) VARNAM(ISQUASH(I)),
     $      (X(I,J),J=JX,MIN(IC,ITOP))
 10   CONTINUE
      JX=JX+5
      IF (JX .GT. IR) RETURN
      WRITE(LuOut,*)
      GOTO 1
 44   FORMAT(T3,A5,T12,F12.6,T22,F12.6,T32,F12.6,T42,
     &       F12.6,T52,F12.6)
#endif 
      RETURN
      END
