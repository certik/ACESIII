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
c get 'NUMVIBRT'
c get TYPE//'SYQT'
c get TYPE//'LABL'
c get TYPE//'DEGN'
c get 'NUMPOINT'
c get 'ENGPOINT'
c get 'INVPSMAT'
c put 'GRADIENT'

      SUBROUTINE SETGRD(NATOM,NIRREP,TYPE,
     &                  LABEL,ISYMIRR,
     &                  SYMGRD,CARTGRD,
     &                  DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 TYPE
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ISYMIRR(3*NATOM)
      DIMENSION SYMGRD(3*NATOM),CARTGRD(3*NATOM)
      double precision dscr(ndscr)

      DIMENSION idegen(100)
      LOGICAL PRINTQ

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
#include "control.com"

      PRINTQ=(IFLAGS(1).GT.10)
      STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5

      NSIZE=3*NATOM
      CALL IGETREC(20,'JOBARC','NUMPOINT',1,NPOINT)
      if (ndscr.lt.npoint) then
         print *, '@SETGRD: Insufficient memory.'
         print *, '         need ',npoint,' doubles'
         print *, '         have ',ndscr,' doubles'
         call aces_exit(1)
      end if
      CALL DGETREC(20,'JOBARC','ENGPOINT',NPOINT,DSCR)
#ifdef _ASSERT
      if (gmtryopt) then
c      o only check the points in irrep 1
         CALL IGETREC(20,'JOBARC','NPTIRREP',1,n)
      else
         n = npoint
      end if
      do i = 1, n
         if (dscr(i).eq.0.d0) then
c         o if any energy is exactly 0, then ACES did not do all points
            print *, '@SETGRD: Assertion failed.'
            print *, '         Energy of point ',i,' is 0. a.u.'
            call aces_exit(1)
         end if
      end do
#endif
      CALL IGETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      CALL IGETREC(20,'JOBARC','INVPSMAT',1,INVOP)
      CALL IGETREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISYMIRR)
      CALL IGETREC(20,'JOBARC',TYPE//'DEGN',NIRREP,IDEGEN)
      IF (PRINTQ) THEN
      CALL DGETREC(20,'JOBARC',TYPE//'LABL',NIRREP,LABEL)
      END IF

      CALL ZERO(SYMGRD,3*NATOM)

      IRREP=1

c   o find first occurance of this irrep
      ILOC=ISRCHEQ(NMODE,ISYMIRR,1,IRREP)
      IF (ILOC.NE.NMODE+1) THEN
         ILAST=ISRCHNE(NMODE,ISYMIRR(ILOC),1,IRREP)
         NVIBSYM=ILAST-1
         NVIBUNQ=NVIBSYM/IDEGEN(IRREP)
         IF (PRINTQ) THEN
            WRITE(6,2000)LABEL(IRREP),IDEGEN(IRREP),NVIBUNQ
2000        FORMAT(T3,' Symmetry : ',A,' Degeneracy : ',I1,
     &             ' Unique symmetry coordinates : ',I3)
         END IF
         if (invop.gt.0) then
            print *, '@SETGRD: Assertion failed.'
            print *, '         Gradients are not implemented for',
     &               ' these displacements.'
            call aces_exit(1)
         end if
         CALL ENER2GRD(NVIBSYM,DSCR,SYMGRD,STPSIZ)
      END IF

c   o transform and write the gradient to JOBARC
      CALL TRNGRD(NATOM,SYMGRD,CARTGRD,DSCR,NDSCR,TYPE,PRINTQ)
      CALL DPUTREC(20,'JOBARC','GRADIENT',NSIZE,CARTGRD)

      RETURN
      END

