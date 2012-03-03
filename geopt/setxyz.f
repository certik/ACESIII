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

      SUBROUTINE SETXYZ(Q,NEWQ,ORIENT,szCOORD,IVAL,IOK,NATOMS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION NEWQ
      CHARACTER*1 szCOORD
      DIMENSION Q(*),NEWQ(*),ORIENT(3,3)
      PARAMETER (ONE=1.0D0)

      IOK=0
      IF (szCOORD.EQ.'X') THEN
         IF (Mod(IVAL,2).EQ.1) THEN
            DO J = 1,NATOMS
               Q(3*J) = NEWQ(3*J-2)
               Q(3*J-2) = NEWQ(3*J-1)
               Q(3*J-1) = NEWQ(3*J)
               ORIENT(1,2)=ONE
               ORIENT(2,3)=ONE
               ORIENT(3,1)=ONE
               IOK=1
            END DO
         END IF
      ELSE
         IF (szCOORD.EQ.'Y') THEN
            IF (Mod(IVAL/2,2).EQ.1) THEN
               DO J = 1,NATOMS
                  Q(3*J) = NEWQ(3*J-1)
                  Q(3*J-1) = NEWQ(3*J-2)
                  Q(3*J-2) = NEWQ(3*J)
                  ORIENT(3,2)=ONE
                  ORIENT(2,1)=ONE
                  ORIENT(1,3)=ONE
                  IOK=1
               END DO
            END IF
         ELSE
            IF (szCOORD.EQ.'Z') THEN
               IF (Mod(IVAL/4,2).EQ.1) THEN
                  CALL XCOPY(3*NATOMS,NEWQ,1,Q,1)
                  ORIENT(1,1)=ONE
                  ORIENT(2,2)=ONE
                  ORIENT(3,3)=ONE
                  IOK=1
               END IF
            END IF
         END IF
      END IF

      RETURN
      END

