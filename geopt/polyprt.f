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
C
      SUBROUTINE POLYPRT(COORD,GRAD,HESS,IORDER,NATOMS,IINTFP,IPRINT,
     &                   LCOORD,LENERG,LGRAD,LHESS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     This subroutine creates a file called POLYRATE and writes to it
C     the Cartesian coordinates, energy, gradient, and hessian, which
C     are read from JOBARC. The units are atomic and the order of atoms
C     is computational.
C
C     TO BE UPDATED :
C     It is
C     called from GEOPT at the time that a vibrational analysis is
C     performed. The main limitation is that it probably will not work
C     when dummy atoms are present.
C
CEND
C
      LOGICAL YESNO,LCOORD,LENERG,LGRAD,LHESS
      DIMENSION COORD(3,NATOMS),GRAD(3,NATOMS),HESS(3*NATOMS,3*NATOMS)
      DIMENSION IORDER(NATOMS)
C
C     LCOORD --- Tells whether to read and write Cartesian coordinates
C     LENERG --- Tells whether to read and write energy
C     LGRAD  --- Tells whether to read and write gradient
C     LHESS  --- Tells whether to read and write hessian
C
C     Make sure this routine has got the correct message about what it is
C     supposed to do.
C
      WRITE(6,3000)
 3000 FORMAT(/,' ----- Entering POLYPRT ----- ',/)
      IF(     LCOORD) WRITE(6,3010)
      IF(.NOT.LCOORD) WRITE(6,3020)
      IF(     LENERG) WRITE(6,3030)
      IF(.NOT.LENERG) WRITE(6,3040)
      IF(     LGRAD ) WRITE(6,3050)
      IF(.NOT.LGRAD ) WRITE(6,3060)
      IF(     LHESS ) WRITE(6,3070)
      IF(.NOT.LHESS ) WRITE(6,3080)
 3010 FORMAT('  @POLYPRT-I, Coordinates will be written. ')
 3020 FORMAT('  @POLYPRT-I, Coordinates will not be written. ')
 3030 FORMAT('  @POLYPRT-I, Energy      will be written. ')
 3040 FORMAT('  @POLYPRT-I, Energy      will not be written. ')
 3050 FORMAT('  @POLYPRT-I, Gradient    will be written. ')
 3060 FORMAT('  @POLYPRT-I, Gradient    will not be written. ')
 3070 FORMAT('  @POLYPRT-I, Hessian     will be written. ')
 3080 FORMAT('  @POLYPRT-I, Hessian     will not be written. ')
C
      IF(LCOORD)THEN
        CALL DGETREC(20,'JOBARC','COORD'   ,3*NATOMS,COORD)
      ENDIF
C
      IF(LENERG)THEN
        CALL DGETREC(20,'JOBARC','TOTENERG',1,ENRG)
      ENDIF
C
      IF(LGRAD)THEN
        CALL DGETREC(20,'JOBARC','GRADIENT',3*NATOMS,GRAD)
      ENDIF
C
      IF(LHESS)THEN
        CALL DGETREC(20,'JOBARC','HESSIANM',9*NATOMS*NATOMS,HESS)
      ENDIF
C
      CALL IGETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IORDER)
C
      IF(IPRINT.GE.1)THEN
C
      WRITE(6,1000)
 1000 FORMAT(/,' Data to be written to POLYRATE file ')
C
      IF(LCOORD)THEN
C
      WRITE(6,1010)
 1010 FORMAT(/,' Cartesian Coordinates ',/)
C
      DO 10 J=1,NATOMS
      WRITE(6,1020) (COORD(I,IORDER(J)),I=1,3)
   10 CONTINUE
 1020 FORMAT(10X,3F20.10)
C
       ENDIF
C
      IF(LENERG)THEN
C
      WRITE(6,1030)
 1030 FORMAT(/,' Total energy ',/)
C
      WRITE(6,1040) ENRG
 1040 FORMAT(F20.12)
C
      ENDIF
C
      IF(LGRAD)THEN
C
      WRITE(6,1050)
 1050 FORMAT(/,' Molecular Gradient (x,y,z for each atom) ',/)
C
      DO 20 J=1,NATOMS
      WRITE(6,1060) (GRAD(I,J),I=1,3)
   20 CONTINUE
 1060 FORMAT(10X,3F20.10)
C
      ENDIF
C
      IF(LHESS)THEN
C
      WRITE(6,1070)
 1070 FORMAT(/,' Molecular Hessian ',/)
C
      DO 40 JOFF=0,3*NATOMS-3,3
      DO 30 I   =1,3*NATOMS
      WRITE(6,1060) (HESS(I,JOFF + J),J=1,3)
   30 CONTINUE
      WRITE(6,1080)
   40 CONTINUE
 1080 FORMAT(/)
C
      ENDIF
C
      ENDIF
C
C     Get the lowest available unit number.
C
      IUNIT=0
      DO 100 I=1,99
C
      IF(IUNIT.EQ.0)THEN
        INQUIRE(UNIT=I,OPENED=YESNO)
          IF(.NOT.YESNO)THEN
            IUNIT = I
            WRITE(6,2000) IUNIT
          ENDIF
      ENDIF
  100 CONTINUE
C
      IF(IUNIT.EQ.0)THEN
      WRITE(6,2010)
      CALL ERREX
      ENDIF
C
C     If POLYRATE already exists, remove it.
C
      INQUIRE(FILE='POLYRATE',EXIST=YESNO)
      IF(YESNO)THEN
      OPEN(UNIT=IUNIT,FILE='POLYRATE',STATUS='OLD',ACCESS='SEQUENTIAL',
     &     FORM='FORMATTED')
      CLOSE(UNIT=IUNIT,STATUS='DELETE')
      ENDIF 
C
C     Create fresh POLYRATE file for current coordinates, energy,
C     gradient, and hessian.
C
      OPEN(UNIT=IUNIT,FILE='POLYRATE',STATUS='NEW',ACCESS='SEQUENTIAL',
     &     FORM='FORMATTED')
C
      WRITE(IUNIT,1080)
C
      IF(LCOORD)THEN
C
      DO 110 J=1,NATOMS
      WRITE(IUNIT,1020) (COORD(I,IORDER(J)),I=1,3)
  110 CONTINUE
C
      WRITE(IUNIT,1080)
C
      ENDIF
C
      IF(LENERG)THEN
C
      WRITE(IUNIT,1040) ENRG
      WRITE(IUNIT,1080)
C
      ENDIF
C
      IF(LGRAD)THEN
C
      DO 120 J=1,NATOMS
      WRITE(IUNIT,1020) (GRAD(I,J),I=1,3)
  120 CONTINUE
C
      WRITE(IUNIT,1080)
C
      ENDIF
C
      IF(LHESS)THEN
C
      DO 140 JOFF=0,3*NATOMS-3,3
      DO 130 I   =1,3*NATOMS
      WRITE(IUNIT,1020) (HESS(I,JOFF + J),J=1,3)
  130 CONTINUE
      WRITE(IUNIT,1080)
  140 CONTINUE
C
      ENDIF
C
      CLOSE(UNIT=IUNIT,STATUS='KEEP')
C
      RETURN
 2000 FORMAT(/,
     & ' @POLYPRT-I, Unit ',I3,' will be used for POLYRATE file. ',/)
 2010 FORMAT(/,' @POLYPRT-F, No available unit for POLYRATE ! ')
      END
