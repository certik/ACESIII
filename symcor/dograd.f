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

c OUTPUT
c double  POINTS(3*NATOM,*) : COORDINATES USED IN FD GRADIENT CALCULATIONS
c integer NPOINT            : NUMBER OF GRADIENT CALCS REQ'D FOR THIS IRREP
c double  DSCR(NDSCR)       : (scr) double scratch

      SUBROUTINE DOGRAD(NATOM,NDIM,SYMQ,COORD,STPSIZ,VMASS,
     &                  POINTS,NPOINT,
     &                  INVOP,PRINTQ,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SYMQ(3*NATOM,NDIM),COORD(3*NATOM),VMASS(NATOM)
      DIMENSION POINTS(*),INVOP(NDIM),DSCR(NDSCR)
      LOGICAL PRINTQ

      LOGICAL ONEGRD
      CHARACTER*5 PHASE

      NSIZE=3*NATOM

      if (ndscr.lt.nsize) then
         print *, '@DOGRAD: Insufficient memory.'
         print *, '         have ',ndscr,' doubles'
         print *, '         need ',nsize,' doubles'
         call aces_exit(1)
      end if

      IOFFP=1
      NPOINT=0

c   o loop over dimensionality of subspace
      DO IDIM=1,NDIM
         ONEGRD=(INVOP(IDIM).GT.0)

c      o scale symmetry coordinate vector
         CALL XDCOPY(NSIZE,SYMQ(1,IDIM),1,DSCR,1)
         CALL XDSCAL(NSIZE,STPSIZ,DSCR,1)

c      o transform from mass-weighted cartesians to pure cartesians
         NDX = 1
         DO IATOM=1,NATOM
            DSCR(NDX+0) = DSCR(NDX+0)*VMASS(IATOM)
            DSCR(NDX+1) = DSCR(NDX+1)*VMASS(IATOM)
            DSCR(NDX+2) = DSCR(NDX+2)*VMASS(IATOM)
            NDX = NDX+3
         END DO

c      o generate positive displacement
         CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,1.d0)
         IOFFP=IOFFP+NSIZE
         NPOINT=NPOINT+1

c      o generate negative displacement (totally symmetric irrep only)
         IF (.NOT.ONEGRD) THEN
            CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,-1.d0)
            IOFFP=IOFFP+NSIZE
            NPOINT=NPOINT+1
         END IF

      END DO

      IF (PRINTQ) THEN
         IOFF=1
         DO IPOINT=1,NPOINT
            IF (.NOT.ONEGRD.AND.MOD(IPOINT,2).EQ.0) THEN
               PHASE='minus'
               ISYCOR=1+(IPOINT-1)/2
            ELSE IF (.NOT.ONEGRD.AND.MOD(IPOINT,2).EQ.1) THEN
               ISYCOR=(IPOINT+1)/2
               PHASE='plus '
            ELSE IF (ONEGRD) THEN
               ISYCOR=IPOINT
               PHASE='plus'
            END IF
            WRITE(6,1001)ISYCOR,PHASE
1001        FORMAT(T3,'Symmetry coordinate : ',i3,' Phase : ',a,
     &             ' Type : Gradient')
            DO IATOM=1,NATOM
               WRITE(6,1002)IATOM,(POINTS(IPOS),IPOS=IOFF,IOFF+2)
1002           FORMAT(T3,I5,3F20.10)
               IOFF=IOFF+3
            END DO
         END DO
      END IF

      RETURN
      END

