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
      SUBROUTINE DUMPCORD(NATOMS,Q,IATNUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mxatms.par"
#include "cbchar.com"
      DIMENSION Q(3*NATOMS),IATNUM(NATOMS)
      ATOB = 0.529177249D0
      WRITE(6,100)
      WRITE(6,101)
      WRITE(6,100)
100   FORMAT(T2,64('-'))
101   FORMAT(T2,'Z-matrix',T13,'Atomic',T31,'C o o r d i n a t e s',
     &       /,
     &       T3,'Symbol',T13,'Number',T30,'X',T45,'Y',T60,'Z')
      IOFF=1
      DO 10 I=1,NATOMS
       LAST=linblnk(ZSYM(I))
       WRITE(6,1000)ZSYM(I)(1:LAST),IATNUM(I),(Q(J)*ATOB,
     &              J=IOFF,IOFF+2)
       IOFF=IOFF+3
10    CONTINUE
      WRITE(6,100)
1000  FORMAT(T6,A,T14,I3,T22,F14.8,T37,F14.8,T52,F14.8)
      RETURN
      END
