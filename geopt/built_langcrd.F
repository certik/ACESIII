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
      SUBROUTINE BULT_LANGCRD(CARTCOORD, BMATRX, ANGL, ICON1, ICON2, 
     &                        ICON3, ICON4, IANGS, NRATMS, TOTREDNCO,
     &                        EPSILON)
C
C Setup the angle bending B-matrix elements for linear bonds. In this case
C we must define two angle bend coordinates.
C
C Assuming that the projection is done on XY-plane. The Z-axis is the 
C perpendicular axis. Given r_i = r_ba and r_j = r_bc then SinG_ij =  
C (x_1*y_2 - y_1*x_2)/R_i*R_j
C
C B(*,a,*,*)(x) =0.0, B(*,*,*,c)(x) = 0.0, 
C B(*,*,b,*) = B(*,a,*,*) - B(*,*,*,c)
C B(*,a,*,*)(y) = -(z_a-z_b)/R_BA, B(*,a,*,*)(z) =  (y_a-y_b)/R_BA
C B(*,*,*,c)(y) =  (z_c-z_b)/R_BC, B(*,*,*,c)(z) = -(y_c-y_b)/R_BC
C B(*,*,b,*)(y) = B(*,a,*,*)(y) + B(*,*,*,c)(y)
C B(*,*,b,*)(z) = B(*,a,*,*)(z) + B(*,*,*,c)(z)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER TOTREDNCO
      DIMENSION CARTCOORD(3*NRATMS), BMATRX(TOTREDNCO, 3*NRATMS),
     &          VECBA(3), VECBC(3), VECAC(3), VECX(3), VECY(3), 
     &          VECZ(3), VECBCX(3), VECBAX(3),  VECN2AC(3),
     &          VECBN2A(3), VECBN2B(3), VECW(3)

C
      DATA VECX /1,0,0/, VECY /0,1,0/, VECZ /0,0,1/
      DATA ZERO /0.0D0/
      DATA MONE /-1/
C
#ifdef _DEBUG_LVL0
          Write(6,*)
          Write(6, "(a)") "Entering the special linear angle block"
#endif
C
      DISTBA = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
      DISTBC = DIST(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2))
 
      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON1 - 2), 
     &         VECBA, 1)
      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECBC, 1)
      CALL VEC(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECAC, 1)     
C
      CALL NORMAL(VECBA, 3)
      CALL NORMAL(VECBC, 3)
      CALL CROSS(VECBA, VECBC, VECN2AC, 0)
      CALL NORMAL(VECN2AC, 3)

#ifdef _DEBUG_LVL0
          Write(6, *)
          Write(6, *) "The normal vector"
          Write(6, "(a, 3F10.4)") "VECN2AC", (VECN2AC(I), I=1,3)
#endif
C
      CALL CROSS(VECBC, VECN2AC, VECW, 1)
      CALL CROSS(VECBC, VECW, VECBN2A, 1)
      CALL CROSS(VECW, VECBA, VECBN2B, 1)
C
#ifdef _DEBUG_LVL0
          Write(6, *)
          Write(6, *) "The basis vctors for special angles"
          Write(6, "(a, 3F10.4)") "VECBN2A", (VECBN2A(I), I=1,3)
          Write(6, "(a, 3F10.4)") "VECBN2B", (VECBN2B(I), I=1,3)
          Write(6,*)
#endif
C
      DO IXYZ = 1, 3
         BMATRX(IANGS, (3*ICON3 - 3) + IXYZ) =
     &                                 VECBN2A(IXYZ)/DISTBC
         BMATRX(IANGS, (3*ICON2 - 3) + IXYZ) =
     &                               - VECBN2A(IXYZ)/DISTBC
     &                               - VECBN2B(IXYZ)/DISTBA
         BMATRX(IANGS, (3*ICON1 - 3) + IXYZ) =
     &                                 VECBN2B(IXYZ)/DISTBA
      ENDDO
C
C Assign the perpendicular bond angle
C

      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
      ANGL = 180.0D0*DINVPI

#ifdef _USE_OLD
      IAXIX = 0
      PRON_X = DABS(XDOT(3, VECX, 1, VECAC, 1))
      PRON_Y = DABS(XDOT(3, VECY, 1, VECAC, 1))
      PRON_Z = DABS(XDOT(3, VECZ, 1, VECAC, 1))
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(A)") "The projections along X, Y and Z axis"
      Write(6, "(3F10.5)") PRON_X, PRON_Y, PRON_Z
#endif
C
      IF (PRON_X .LT. EPSILON) THEN
C
C The molecule is in YZ plane, X-axis is the perpendicular axis for 
C one of the two angles
C 
         IF (PRON_Y .LT. EPSILON) THEN
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 1 
            ELSE
               IAXIX = 2
            ENDIF
         ENDIF
C         
         IF (PRON_Z .LT. EPSILON) THEN
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 1
            ELSE
               IAXIX = 3
            ENDIF
         ENDIF
C
      ELSE IF (PRON_Y .LT. EPSILON) THEN
C
C The molecule is in XZ plane, Y-axis  is the perpendicular axis for 
C one of the two angles
C
         IF (PRON_X .LT. EPSILON) THEN
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 2
            ELSE
               IAXIX = 1
            ENDIF
         ENDIF
C         
         IF (PRON_Z .LT. EPSILON) THEN
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 2
            ELSE
               IAXIX = 3
            ENDIF
         ENDIF
C
      ELSE IF (PRON_Z .LT. EPSILON) THEN
C
C The molecule is in XY plane Z-axis is the perpendicular axis for 
C one of the two angles
C
         IF (PRON_X .LT. EPSILON) THEN
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 3
            ELSE
               IAXIX = 1
            ENDIF
         ENDIF
C         
         IF (PRON_Y .LT. EPSILON) THEN
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 3
            ELSE
               IAXIX = 2 
            ENDIF
         ENDIF
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a,a,i)") "The Iaxix after the projection logic",
     &                   " level 1:", IAXIX
      Write(6,*)
#endif
C
      ELSE
C
C Molecule is not in any of the above planes. We have to look at the 
C size of the projections on XY, YZ, XZ planes to decide axises.
C  
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a,3F10.5)") "The AC vector",(VECAC(I), I=1,3)
      Write(6,*)
#endif
C     
         IF (DABS(VECAC(1)) .LT. DABS(VECAC(2)) .AND. DABS(VECAC(1))
     &       .LT. DABS(VECAC(3))) THEN            
C 
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 1
            ELSE
               IF (ABS(VECAC(2)) .LT. ABS(VECAC(3))) THEN
                  IAXIX = 2
               ENDIF
            ENDIF
C
         ELSE IF (ABS(VECAC(2)) .LT. ABS(VECAC(1)) .AND. ABS(VECAC(2))
     &            .LT. ABS(VECAC(3))) THEN            
C
            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 2
            ELSE
               IF (ABS(VECAC(1)) .LT. ABS(VECAC(3))) THEN
                  IAXIX = 1
               ENDIF
            ENDIF
C
         ELSE IF (ABS(VECAC(3)) .LT. ABS(VECAC(1)) .AND. ABS(VECAC(3)) 
     &            .LT. ABS(VECAC(2))) THEN            

            IF (ICON4 .EQ. MONE) THEN
               IAXIX = 3
            ELSE
               IF (ABS(VECAC(1)) .LT. ABS(VECAC(2))) THEN
                  IAXIX = 1
               ENDIF
            ENDIF
C
         ENDIF
C
      ENDIF
C
C Get the projection vectors on appropriate planes, That means setting the 
C the component of BA and BC vectors along the projection axis to zero!
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a,i)") "The Iaxix after the projection logic:", IAXIX
      Write(6,*)
#endif     
      VECBA(IAXIX) = ZERO
      VECBC(IAXIX) = ZERO
      ANGL         = ZERO 
      DISTBA = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
      DISTBC = DIST(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2))
C     
      IF (IAXIX .EQ. 1) THEN
C
         CALL CROSS(VECX, VECBA, VECBAX, 0) 
         CALL CROSS(VECX, VECBC, VECBCX, 0) 
         BMATRX(IANGS, 3*ICON1 + 2) = -VECBAX(3)/DISTBA
         BMATRX(IANGS, 3*ICON3 + 2) =  VECBCX(3)/DISTBC
         BMATRX(IANGS, 3*ICON2 + 2) = -VECBAX(3)/DISTBA
     &                                +VECBCX(3)/DISTBC 
         BMATRX(IANGS, 3*ICON1 + 3) =  VECBAX(2)/DISTBA
         BMATRX(IANGS, 3*ICON3 + 3) = -VECBCX(2)/DISTBC
         BMATRX(IANGS, 3*ICON2 + 3) =  VECBAX(2)/DISTBA
     &                                -VECBCX(2)/DISTBC
C
      ELSE IF (IAXIX .EQ. 2) THEN
C
         CALL CROSS(VECY, VECBA, VECBAX, 0) 
         CALL CROSS(VECY, VECBC, VECBCX, 0) 
         BMATRX(IANGS, 3*ICON1 + 1) =  VECBAX(3)/DISTBA
         BMATRX(IANGS, 3*ICON3 + 1) = -VECBCX(3)/DISTBC
         BMATRX(IANGS, 3*ICON2 + 1) =  VECBAX(3)/DISTBA  
     &                                -VECBCX(3)/DISTBC
         BMATRX(IANGS, 3*ICON1 + 3) = -VECBAX(1)/DISTBA
         BMATRX(IANGS, 3*ICON3 + 3) =  VECBCX(1)/DISTBC
         BMATRX(IANGS, 3*ICON2 + 3) = -VECBAX(1)/DISTBA  
     &                                +VECBCX(1)/DISTBC
C
      ELSE IF (IAXIX .EQ. 3) THEN
C
         CALL CROSS(VECZ, VECBA, VECBAX, 0) 
         CALL CROSS(VECZ, VECBC, VECBCX, 0) 
         BMATRX(IANGS, 3*ICON1 + 1) =  VECBAX(2)/DISTBA
         BMATRX(IANGS, 3*ICON3 + 1) = -VECBCX(2)/DISTBC 
         BMATRX(IANGS, 3*ICON2 + 1) =  VECBAX(2)/DISTBA
     &                                -VECBCX(2)/DISTBC 
         BMATRX(IANGS, 3*ICON1 + 2) = -VECBAX(1)/DISTBA
         BMATRX(IANGS, 3*ICON3 + 2) =  VECBCX(1)/DISTBC
         BMATRX(IANGS, 3*ICON2 + 2) = -VECBAX(1)/DISTBA  
     &                                +VECBCX(1)/DISTBC
C
      ENDIF
#endif
C
      RETURN
      END
         
