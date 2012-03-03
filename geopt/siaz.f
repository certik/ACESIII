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
      SUBROUTINE SIAZ(VEC,ROT,IAXIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(3),ROT(3,3),RM(3,3)
C 
C DETERMINE LENGTH OF VECTOR.
C
      DIS=DSQRT(xdot(3,VEC,1,VEC,1))
C
C RETURNS ROTATION MATRIX NEEDED TO PLACE VECTOR
C   VEC ALONG THE X (IAXIS=1), Y (IAXIS=2) OR Z (IAXIS=3)
C   AXIS.  TAIL OF VECTOR IMPLICITLY ASSUMED TO BE AT 
C   THE ORIGIN.
C
      IF(IAXIS.EQ.3)THEN
        ARG1=VEC(3)/DIS
        IF(DABS(ARG1).GT.1.D0)ARG1=ARG1/DABS(ARG1)
        AZIM=-1.D0*DACOS(ARG1)
        DISP=DSQRT(xdot(2,VEC,1,VEC,1))
        IF(ABS(DISP).LT.1.D-10)THEN
         CALL ROTM(3,0.0D0,1,ROT)
         RETURN
        ENDIF
        ARG2=VEC(1)/DISP
        IF(DABS(ARG2).GT.1.D0)ARG2=ARG2/DABS(ARG2)
        IF(DABS(DABS(ARG2)-1.D0).LT.1.D-10)ARG2=SIGN(1.D0,ARG2)
        ANG1=DACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(2))
        CALL ROTM(3,ANG1,0,RM)
        CALL ROTM(2,AZIM,0,ROT)
        CALL MODMATMUL(ROT,RM,ROT,3,3,3,3,3,3)
      ELSEIF(IAXIS.EQ.1)THEN
        ARG1=VEC(1)/DIS
        IF(DABS(ARG1).GT.1.D0)ARG1=ARG1/DABS(ARG1)
        AZIM=-1.D0*DACOS(ARG1)
        DISP=DSQRT(xdot(2,VEC(2),1,VEC(2),1))
        IF(ABS(DISP).LT.1.D-10)THEN
         CALL ROTM(3,0.0D0,1,ROT)
         RETURN
        ENDIF
        ARG2=VEC(2)/DISP
        IF(DABS(ARG2).GT.1.D0)ARG2=ARG2/DABS(ARG2)
        IF(DABS(DABS(ARG2)-1.D0).LT.1.D-10)ARG2=SIGN(1.D0,ARG2)
        ANG1=DACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(3))
        CALL ROTM(1,ANG1,0,RM)
        CALL ROTM(3,AZIM,0,ROT)
        CALL MODMATMUL(ROT,RM,ROT,3,3,3,3,3,3)
      ELSEIF(IAXIS.EQ.2)THEN
        ARG1=VEC(2)/DIS
        IF(DABS(ARG1).GT.1.D0)ARG1=ARG1/DABS(ARG1)
        AZIM=-1.D0*DACOS(ARG1)
        DISP=DSQRT(DIS**2-VEC(2)**2)
        IF(ABS(DISP).LT.1.D-10)THEN
         CALL ROTM(3,0.0D0,1,ROT)
         RETURN
        ENDIF
        ARG2=VEC(3)/DISP
        IF(DABS(ARG2).GT.1.D0)ARG2=ARG2/DABS(ARG2)
        IF(DABS(DABS(ARG2)-1.D0).LT.1.D-10)ARG2=SIGN(1.D0,ARG2)
        ANG1=DACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(1))
        CALL ROTM(2,ANG1,0,RM)
        CALL ROTM(1,AZIM,0,ROT)
        CALL MODMATMUL(ROT,RM,ROT,3,3,3,3,3,3)
       ENDIF
       RETURN
       END
