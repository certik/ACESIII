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
c double  DSCR(*)
c integer ISCR(*)
c char*8  LABEL(NIRREP+NIRREPC)

c RECORDS
c get DOIT//'SYMQ'
c get DOIT//'SYQT'
c get DOIT//'LABL'
c get 'NUMVIBRT'
c get 'LINEAR  '
c get 'SBGRPSYM'
c get 'COMPLABL'
c get DOIT//'PTGP'
c get 'COMPPTGP'

      SUBROUTINE PRTCOORD(NATOM,NIRREP,NIRREPC,DSCR,ISCR,
     &                    LABEL,DOIT,BFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION DSCR(*),ISCR(*)
      CHARACTER*8 LABEL(NIRREP+NIRREPC)
      CHARACTER*4 DOIT
      LOGICAL BFLAG

      CHARACTER*4 FPGRP,CPGRP

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "control.com"

      NSIZE=3*NATOM

      CALL DGETREC(20,'JOBARC',DOIT//'SYMQ',NSIZE*NSIZE,DSCR)
      CALL IGETREC(20,'JOBARC',DOIT//'SYQT',NSIZE,ISCR)
      CALL DGETREC(20,'JOBARC',DOIT//'LABL',NIRREP,LABEL)
      CALL IGETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      CALL IGETREC(20,'JOBARC','LINEAR  ',1,ILINEAR)
      IF (BFLAG) THEN
      CALL IGETREC(20,'JOBARC','SBGRPSYM',NSIZE,ISCR(NSIZE+1))
      CALL DGETREC(20,'JOBARC','COMPLABL',NIRREPC,LABEL(NIRREP+1))
      CALL GETCREC(20,'JOBARC',DOIT//'PTGP',4,FPGRP)
      CALL GETCREC(20,'JOBARC','COMPPTGP',  4,CPGRP)
      END IF

      IOFF=1
      WRITE(6,1000)
1000  FORMAT(T3,' *** TABLE OF SYMMETRY COORDINATES *** ')
      IF (BFLAG) THEN
      WRITE(6,1010)
1010  FORMAT(T3,'                                       ')
      IF (ROTPROJ) THEN
         WRITE(6,1020)
1020     FORMAT(T3,' Eckart rotations and translations projected ',
     &             'from vibrational coordinates')
      ELSE
         WRITE(6,1021)
1021     FORMAT(T3,' Translations projected from vibrational ',
     &             'coordinates')
      END IF
      WRITE(6,1030)FPGRP,CPGRP
1030  FORMAT(T3,' Point group used : ',A,' ; Abelian subgroup : ',A)
      WRITE(6,1031)DOIT
1031  FORMAT(T3,'Coordinate type : ',A)
c     END IF (BFLAG)
      END IF
      DO IMODE=1,NMODE
         IF (BFLAG) THEN
         WRITE(6,1003)IMODE,LABEL(ISCR(IMODE)),
     &                LABEL(NIRREP+ISCR(NSIZE+IMODE))
1003     FORMAT(T3,' Mode : ',I3,' Symmetry : ',A,' Subgroup symm : ',A)
         ELSE
         WRITE(6,1001)IMODE,LABEL(ISCR(IMODE))
1001     FORMAT(T3,' Mode : ',I3,' Symmetry : ',A)
         END IF
         WRITE(6,1004)
1004     FORMAT(T11,'ATOM',T26,'X',T40,'Y',T54,'Z')
         DO IATOM=1,NATOM
            WRITE(6,1005)IATOM,(DSCR(J),J=IOFF,IOFF+2)
            IOFF=IOFF+3
1005        FORMAT(T11,I3,T19,F13.10,T33,F13.10,T47,F13.10)
         END DO
         WRITE(6,*)
      END DO

      IF (ROTPROJ) THEN
         WRITE(6,2000)
2000     FORMAT(T3,'@PRTCOORD: Rotational coordinates (x,y,z)')
         DO IXYZ=1,3-ILINEAR
            WRITE(6,1004)
            DO IATOM=1,NATOM
               WRITE(6,1005)IATOM,(DSCR(J),J=IOFF,IOFF+2)
               IOFF=IOFF+3
            END DO
         END DO
      END IF

      WRITE(6,3000)
3000  FORMAT(T3,'@PRTCOORD: Translational coordinates (x,y,z)')
      DO IXYZ=1,3
         WRITE(6,1004)
         DO IATOM=1,NATOM
            WRITE(6,1005)IATOM,(DSCR(J),J=IOFF,IOFF+2)
            IOFF=IOFF+3
         END DO
      END DO

      RETURN
      END

