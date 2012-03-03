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
c get TYPE//'SYMQ'
c get 'ATOMMASS'

      SUBROUTINE TRNGRD(NATOM,SYMGRD,CARTGRD,DSCR,NDSCR,TYPE,PRINTQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION DSCR(NDSCR)
      DOUBLE PRECISION SYMGRD(3*NATOM),CARTGRD(3*NATOM)
      CHARACTER*4 TYPE
      LOGICAL PRINTQ

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      NSIZE=3*NATOM

      IF (NDSCR.LT.NSIZE*NSIZE) THEN
         print *, '@TRNGRD: Insufficient memory to load ',TYPE,'SYMQ'
         print *, '         have ',NDSCR,' doubles'
         print *, '         need ',NSIZE*NSIZE,' doubles'
         call aces_exit(1)
      END IF

c   o transform to mass-weighted cartesian coordinates
      CALL DGETREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE,DSCR)
      CALL XGEMM('N','N',NSIZE,1,NSIZE,
     &           1.d0,DSCR,   NSIZE,
     &                SYMGRD, NSIZE,
     &           0.d0,CARTGRD,NSIZE)

c   o remove mass weighting
      CALL DGETREC(20,'JOBARC','ATOMMASS',NATOM,DSCR)
      IOFF=1
      DO IATOM=1,NATOM
         X=SQRT(DSCR(IATOM))
         CARTGRD(IOFF+0)=X*CARTGRD(IOFF+0)
         CARTGRD(IOFF+1)=X*CARTGRD(IOFF+1)
         CARTGRD(IOFF+2)=X*CARTGRD(IOFF+2)
         IOFF=IOFF+3
      END DO

      IF (PRINTQ) THEN
         write(6,*)' Gradient before transformation:'
         write(6,'((3f20.10))')(symgrd(i),i=1,nsize)
         write(6,*)' Numerical Cartesian gradient:'
         write(6,'((3f20.10))')(cartgrd(i),i=1,nsize)
      END IF

      RETURN
      END

