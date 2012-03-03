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
      SUBROUTINE GETICFCM(A,NIC,IICFCM)
C
C ROUTINE FETCHES INTERNAL OR CARTESIAN COORDINATE HESSIAN FROM FCMINT.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL YESNO
#include "fnamelen.par"
      CHARACTER*(fnamelen) FNAME
      DIMENSION A(NIC,NIC)
      IICFCM=0
      CALL GFNAME('FCMINT  ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
      IF(.NOT.YESNO)RETURN
c      OPEN(UNIT=50,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
      OPEN(UNIT=50,FILE=FNAME(1:ILENGTH),FORM='FORMATTED',
     &     STATUS='OLD')
c      REWIND(50)
      READ(50,*)((A(I,J),J=1,NIC),I=1,NIC)
      CLOSE(UNIT=50,STATUS='KEEP')
      WRITE(6,100)
100   FORMAT(T3,'@GETICFCM-I, Hessian read from FCMINT file.')
      IICFCM=1
      RETURN
      END
