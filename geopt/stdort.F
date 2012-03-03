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
      SUBROUTINE STDORT(VEC,SCRATCH,LEN,ICNTL)
C
C SUBROUTINE ROTATES MOLECULE TO THE PRINCIPAL AXIS ORIENTATION FOR
C  THE FULL MOLECULAR POINT GROUP.  ICNTL=0 ROTATES THE MOLECULE FROM
C  THE ACES2 ORIENTATION TO THE PAO, ANY OTHER VALUE ROTATES FROM THE
C  PAO TO ACES2 ORIENTATION.  LEN IS THE LENGTH OF THE VECTOR PASSED
C  IN.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "io_units.par"
      Common /Orient/ Orient(3,3)
      DIMENSION VEC(LEN),SCRATCH(LEN)
C
#ifdef _DEBUG_LVLM1
      Print*, "The orient in stdort", Orient 
#endif
c
      NTIMES=LEN/3
      IOK=0
      IF(MOD(LEN,3).NE.0)THEN
       WRITE(LUOUT,133)LEN
 133   FORMAT(T3,'@STDORT-F, Vector with strange length [',i6,'] passed'
     &' to STDORT.')
       Call ErrEx
      ENDIF
      DO 30 I=1,3
       DO 31 J=1,3
       IF(ABS(ORIENT(I,J)).GT.1.D-3)IOK=IOK+1
 31    CONTINUE
 30   CONTINUE
      IF(IOK.EQ.0)THEN
       WRITE(6,5432)
 5432 FORMAT(T3,'@STDORT-W, Orientation matrix is empty.  Unit ',
     &'matrix used.')
      DO 40 I=1,3
       ORIENT(I,I)=1.D0
 40   CONTINUE
      ENDIF
      CALL ZERO (SCRATCH,LEN)
      CALL VADD (SCRATCH,SCRATCH,VEC,LEN,1.D0)
      CALL ZERO (VEC,LEN)
      IF(ICNTL.EQ.0)CALL MTRAN2(ORIENT,3)
      CALL MODMATMUL(VEC,ORIENT,SCRATCH,3,3,NTIMES,3,3,NTIMES)
      IF(ICNTL.EQ.0)CALL MTRAN2(ORIENT,3)
      IF(IOK.EQ.0)CALL ZERO(ORIENT,9)
      RETURN
      END
