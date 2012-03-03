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
c double  GRDXYZ(3,NATOM)
c double  GRDINT(3,NATOM)
c double  DIPXYZ(3)
c double  POLXYZ(3,3)
c integer IMAP(NATOM)
c double  DSCR(NDSCR)

c RECORDS
c get 'COORD'
c get 'ATOMMASS'
c get 'MAP2ZMAT'
c get 'NREALATM'
c get 'GRADIENT'
c get DOIT//'SYMQ'

      SUBROUTINE GETGRD(NATOM,DOIT,REFGEOM,
     &                  GRDXYZ,GRDINT,DIPXYZ,POLXYZ,
     &                  IMAP,DSCR,NDSCR,PRINTQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 DOIT
      DIMENSION REFGEOM(3,NATOM)
      DIMENSION GRDXYZ(3,NATOM),GRDINT(3,NATOM)
      DIMENSION DIPXYZ(3),POLXYZ(3,3)
      DIMENSION IMAP(NATOM)
      double precision dscr(ndscr)
      LOGICAL PRINTQ

      DIMENSION DIPOL(3)
      double precision WMAT(3,3), DA3x3(3,3)
      LOGICAL POLAR_EXIST, DIPOL_EXIST

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      NSIZE = 3*NATOM
      lTMP  = 1+NSIZE

      if (ndscr.lt.nsize*(1+nsize)) then
         print *, '@GETGRD: Insufficient memory.'
         print *, '         have ',ndscr,' doubles'
         print *, '         need ',nsize*(1+nsize),' doubles'
         call aces_exit(1)
      end if

c The scratch usage works like this:
c     COORD(1:NSIZE) -> DSCR(1)
c     ATOMMASS(1:NATOM) -> DSCR(1+NSIZE)
c     [ COORD + ATOMMASS forms WMAT ]
c     [ WMAT + GRADIENT forms GRDSCR(1:NSIZE) in DSCR(1) ]
c     [ GRDSCR + ATOMMASS forms GRDSCR2(1:NSIZE) in DSCR(1) ]
c     DOIT//'SYMQ'(1:NSIZE,1:NSIZE) -> DSCR(1+NSIZE)

c   o get the rotation matrix that maps the computational geometry into the
c     original reference geometry from the finite difference grid
      CALL DGETREC(20,'JOBARC','COORD',NSIZE,DSCR)
      CALL DGETREC(20,'JOBARC','ATOMMASS',NATOM,DSCR(lTMP))
      CALL Q2QPRIME(DSCR,REFGEOM,DSCR(lTMP),DSCR(lTMP+NATOM),WMAT,
     &              NATOM)

C READ GRADIENT (IN VMOL ORDER) AND PUT ATOMS IN CORRECT POSITIONS
C IN GRD VECTOR. CHANGE TO READ DIRECTLY FROM JOBARC 07/2000, Ajith Perera.

C See notes in vmol2ja for ordering reordering issues. Both Vmol and Seward
C generates symmetry redundent atoms in the same order and the MAP2ZMAT
C is same for both cases. See vmol2ja (v2j.f) for detailed notes.

c   o pick up VMOL->ZMAT mapping vector
      CALL IGETREC(20,'JOBARC','MAP2ZMAT',NATOM,IMAP)
      CALL IGETREC(20,'JOBARC','NREALATM',1,NRATOM)
      CALL DGETREC(20,'JOBARC','GRADIENT',3*NRATOM,GRDINT)

      CALL ZERO(GRDXYZ,NSIZE)
      DO IATMVML=1,NATOM
         IATMZMAT=IMAP(IATMVML)
         IF (IATMZMAT.NE.0) THEN
            GRDXYZ(1,IATMZMAT)=GRDINT(1,IATMVML)
            GRDXYZ(2,IATMZMAT)=GRDINT(2,IATMVML)
            GRDXYZ(3,IATMZMAT)=GRDINT(3,IATMVML)
         END IF
      END DO

c   o transform the gradient from int orient to ext orient
      CALL XGEMM('N','N',3,NATOM,3,
     &           1.d0,WMAT,  3,
     &                GRDXYZ,3,
     &           0.d0,DSCR,  3)
      CALL XDCOPY(NSIZE,DSCR,1,GRDXYZ,1)

c   o mass-weight the gradient
      IOFF=1
      DO IATOM=0,NATOM-1
         IF (DSCR(lTMP+IATOM).NE.0.d0) THEN
            X=1.d0/SQRT(DSCR(lTMP+IATOM))
            DSCR(IOFF+0)=X*DSCR(IOFF+0)
            DSCR(IOFF+1)=X*DSCR(IOFF+1)
            DSCR(IOFF+2)=X*DSCR(IOFF+2)
         END IF
         IOFF=IOFF+3
      END DO

c   o transform the gradient to symmetry adapted coordinates
      CALL DGETREC(20,'JOBARC',DOIT//'SYMQ',NSIZE*NSIZE,DSCR(lTMP)
     &           )
      CALL XGEMM('T','N',NSIZE,1,NSIZE,
     &           1.d0,DSCR(lTMP),NSIZE,
     &                DSCR,   NSIZE,
     &           0.d0,GRDINT, NSIZE)

C NOW DEAL WITH THE DIPOLE MOMENT. Note that for the moment
C We can't do IR intensities with ALASKA integrals. We need to
C find a best possible way to read dipole integrals from ONE_INT
C file and the CC density matrix from JOBARC (or lists) to create
C the dipole moment vector. This is a part of the job that was
C assigned to Carlos Taylor. 07/2000, Ajith Perera

      INQUIRE(FILE='DIPOL',EXIST=DIPOL_EXIST)
      IF (DIPOL_EXIST) THEN
         OPEN(UNIT=10,FILE='DIPOL',FORM='FORMATTED',STATUS='OLD')
         READ(10,'(3F20.10)')DIPOL
         CALL XGEMM('N','N',3,1,3,1.d0,WMAT,3,DIPOL,3,0.d0,DIPXYZ,3)
         CLOSE(UNIT=10,STATUS='KEEP')
      END IF

C Now deal with the polarizability, 07/98 John and Ajith

      INQUIRE(FILE='POLAR',EXIST=POLAR_EXIST)
      IF (POLAR_EXIST) THEN
         OPEN(UNIT=10,FILE='POLAR',FORM='FORMATTED',STATUS='OLD')
         DO I = 1, 3
            READ(10,'(3F20.10)') (POLXYZ(J,I),J=1,3)
         END DO
         CLOSE(UNIT=10,STATUS='KEEP')
      END IF

c   o transform to the reference coordinate system
      CALL XGEMM('N','N',3,3,3,1.d0,WMAT,3,POLXYZ,3,0.d0,DA3x3,3)
      CALL XGEMM('N','T',3,3,3,1.d0,DA3x3,3,WMAT,3,0.d0,POLXYZ,3)

      IF (PRINTQ) THEN
         WRITE(6,*)' Dipole moment '
         WRITE(6,'(3F20.10)')(DIPXYZ(I),I=1,3)
         WRITE(6,*)' Gradient vector in cartesian coordinates '
         DO IATOM=1,NATOM
            WRITE(6,'(I5,3F20.10)')IATOM,(GRDXYZ(I,IATOM),I=1,3)
         END DO
         WRITE(6,*)' Gradient vector in internal coordinates '
         DO IATOM=1,NATOM
            WRITE(6,'(I5,3F20.10)')IATOM,(GRDINT(I,IATOM),I=1,3)
         END DO
         WRITE(6,*)' Polarizability '
         DO I=1, 3
            WRITE(6,'(3F20.10)')(POLXYZ(I,J),J=1,3)
         END DO
      END IF

      RETURN
      END

