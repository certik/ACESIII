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
      SUBROUTINE BULT_DIHANGCRD(CARTCOORD, BMATRX, DANG, ICON1, 
     &                          ICON2, ICON3, ICON4, IDIHS, 
     &                          TOTREDNCO, NRATMS)
C
C Setup the dihedral angle B-matrix elements. 
C It is some what tedious to write down the elements of the B-matrix for
C torsion. Please refer to Wilson, Decius and Cross page 61, 
C Eqns. 21,22,23,24.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VECAB , VECBC , VECBA , DEGPI , VECCD
      INTEGER TOTREDNCO
      DIMENSION CARTCOORD(3*NRATMS), BMATRX(TOTREDNCO, 3*NRATMS),
     &  VECCB(3), VECAB(3), VECBC(3), VECCD(3), VECABC(3), VECBCD(3)
C
      DATA DEGPI /180.0D0/
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
C
#ifdef _DEBUG_LVLM1
      WRITE(6,"(4I4)") icon1, icon2, icon3, icon4
#endif
      CALL VEC(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2), 
     &         VECAB, 1)
      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECBC, 1)
      CALL VEC(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2), 
     &         VECCB, 1)
      CALL VEC(CARTCOORD(3*ICON4 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECCD, 1)
C
      DISTAB = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
      DISTBC = DIST(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2))
      DISTCD = DIST(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON4 - 2))
C
C First evaluate the dihedral angle. This is calculated as the 
C angle between the two unit vectors that are perpendicular to 
C the ABC and BCD planes for A-B-C-D pattern. 
C
      CALL CROSS(VECAB, VECBC, VECABC, 1)
C
C We need the vector CD, obtain that from vec DC (note the misnomar)
C
      CALL DSCAL(3, -1.0D0, VECCD, 1)
      CALL CROSS(VECBC, VECCD, VECBCD, 1)
      DANG   = (ANGLE(VECABC, VECBCD, 3))*DINVPI
      
C
C Retain the vector DC for continuation. 
C
      CALL DSCAL(3, -1.0D0, VECCD, 1)
C
C Now re-evalute the unnormalied vectors along the the perpendicular
C directions to the planes ABC and BCD. 
C
      CALL CROSS(VECAB, VECBC, VECABC, 0)
      CALL CROSS(VECCD, VECCB, VECBCD, 0)
C
      DANG1  = ANGLE(VECAB, VECCB, 3)*DINVPI
      DANG2  = ANGLE(VECBC, VECCD, 3)*DINVPI 
C
      SINABC = DSIN(DANG1)
      COSABC = DCOS(DANG1)
      SINBCD = DSIN(DANG2)
      COSBCD = DCOS(DANG2)
C
      DABC = DDOT(3, VECABC, 1, VECABC, 1)
      DBCD = DDOT(3, VECBCD, 1, VECBCD, 1)
C
C--- Obtainig the sign of the dihedral angle:
C
      SENSE = DDOT(3, VECAB, 1, VECBCD, 1)

      IF (SENSE .GT. 0.0D0) DANG = -DANG
C
#ifdef _DEBUG_LVLM2
      Print*, "The vectors making the dihedral angle"
      Print*,  (CARTCOORD(3*ICON1 - 2 + I), I=0, 2)
      Print*,  (CARTCOORD(3*ICON2 - 2 + I), I=0, 2)
      Print*,  (CARTCOORD(3*ICON3 - 2 + I), I=0, 2)
      Print*,  (CARTCOORD(3*ICON4 - 2 + I), I=0, 2)
      Print*, "The vec. BC: ", (VECBC(I), I=1, 3)
      Print*, "The vec. AB: ", (VECAB(I), I=1, 3)
      Print*, "The vec. CB: ", (VECCB(I), I=1, 3)
      Print*, "The vec. CD: ", (-VECCD(I), I=1, 3)
      Print*, "The vec. ABC:", (VECABC(I), I=1,3)
      Print*, "The vec. BCD:", (VECBCD(I), I=1,3)
      Write(6,"(a,F10.5)") "The dihedral angle", DANG/DINVPI
#endif

      IF (DABC .GT. 1.0D-16 .AND. DBCD .GT. 1.0D-16) THEN

      DO 10 IXYZ = 1, 3 
C
         BMATRX(IDIHS, (3*ICON1 - 3) + IXYZ) = 
     &                     -(VECABC(IXYZ)/(DISTAB*SINABC*SINABC))
                                            
         BMATRX(IDIHS, (3*ICON2 - 3) + IXYZ) = 
     &                  ((DISTBC - DISTAB*COSABC)*VECABC(IXYZ)/
     &                        (DISTAB*DISTBC*SINABC*SINABC)) + 
     &                            ((COSBCD*VECBCD(IXYZ))/
     &                            (DISTBC*SINBCD*SINBCD)) 
C
         BMATRX(IDIHS, (3*ICON3 - 3) + IXYZ) = 
     &                     ((DISTBC - DISTCD*COSBCD)*VECBCD(IXYZ)/
     &                         (DISTCD*DISTBC*SINBCD*SINBCD)) + 
     &                                 ((COSABC*VECABC(IXYZ))/
     &                                 (DISTBC*SINABC*SINABC)) 
C
         BMATRX(IDIHS, (3*ICON4 - 3) + IXYZ) = 
     &                     -(VECBCD(IXYZ)/(DISTCD*SINBCD*SINBCD))
         
 10   CONTINUE
      ENDIF 
C
      RETURN
      END


      
