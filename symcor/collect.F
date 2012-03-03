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

c RECORDS
c get TYPE//'SYMQ'
c get TYPE//'SYQT'
c get TYPE//'NSYQ'
c get TYPE//'LABL'
c get 'LINEAR'
c get ROTREC(IXYZ)
c get TRAREC(IXYZ)
c put 'NVIBSYM'
c put TYPE//'SYMQ'
c put TYPE//'SYQT'
c put 'NUMVIBRT'

      SUBROUTINE COLLECT(NATOM,NIRREP,SCRIN,SCROUT,
     &                   ISCRIN,ISCROUT,LABEL,TYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SCRIN(*),SCROUT(*),ISCRIN(*),ISCROUT(*)
      CHARACTER*8 LABEL(NIRREP)
      CHARACTER*4 TYPE

      CHARACTER*8 ROTREC(3),TRAREC(3)
      DIMENSION NVIBSYM(20)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "control.com"

      DATA TOL/1.D-8/
      DATA ROTREC /'ROTVECX ','ROTVECY ','ROTVECZ '/
      DATA TRAREC /'TRAVECX ','TRAVECY ','TRAVECZ '/

      IPOSOUT=0
      NSIZE=3*NATOM

c   o read in symmetry adapted coordinates and their symmetries
      CALL DGETREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE,SCRIN)
      CALL IGETREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISCRIN)
      CALL IGETREC(20,'JOBARC',TYPE//'NSYQ',1,NCOORD)
      CALL DGETREC(20,'JOBARC',TYPE//'LABL',NIRREP,LABEL)
      CALL IGETREC(20,'JOBARC','LINEAR  ',1,ILINEAR)

c   o loop over coordinates and remove zero entries
      CALL IZERO(NVIBSYM,NIRREP)
      DO IPOSIN=1,NCOORD
         LOCIN=1+(IPOSIN-1)*NSIZE
         X=DNRM2(NSIZE,SCRIN(LOCIN),1)
         IF (X.GT.TOL) THEN
            IPOSOUT=IPOSOUT+1
            LOCOUT=1+(IPOSOUT-1)*NSIZE
            CALL XDCOPY(NSIZE,SCRIN(LOCIN),1,SCROUT(LOCOUT),1)
            ISCROUT(IPOSOUT)=ISCRIN(IPOSIN)
            NVIBSYM(ISCROUT(IPOSOUT))=NVIBSYM(ISCROUT(IPOSOUT))+1
         END IF
      END DO

      NMODE=IPOSOUT

      ILOC=IPOSOUT*NSIZE+1

c   o put rotations and translations in rest of matrix
      IF (ROTPROJ) THEN
         ITOP=3-ILINEAR
         DO IXYZ=1,ITOP
            CALL DGETREC(20,'JOBARC',ROTREC(IXYZ),NSIZE,
     &                  SCROUT(ILOC))
            ILOC=ILOC+NSIZE
         END DO
      END IF
      DO IXYZ=1,3
         CALL DGETREC(20,'JOBARC',TRAREC(IXYZ),NSIZE,SCROUT(ILOC))
         ILOC=ILOC+NSIZE
      END DO

      CALL IPUTREC(20,'JOBARC','NVIBSYM ',NIRREP,NVIBSYM)
      CALL DPUTREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NCOORD,SCROUT)
      CALL IPUTREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISCROUT)
      CALL IPUTREC(20,'JOBARC','NUMVIBRT',1,NMODE)

      RETURN
      END

