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
c integer IMORE
c double  ENGPT(9*NATOM*NATOM)
c double  GRDPT(3*NATOM,9*NATOM*NATOM)
c double  DIPPT(3,9*NATOM*NATOM)
c double  POLPT(9,9*NATOM*NATOM)
c integer IPTTYPE(9*NATOM*NATOM)
c double  DSCR(NDSCR)

c RECORDS
c get 'NUMPOINT'
c get 'FDCALCTP'
c get 'FDCOORDS'
c put 'NEXTGEOM'
c put 'FDCALCTP'
c get 'ENGPOINT'
c get 'TOTENERG'
c get 'TOTENER2'
c put 'ENGPOINT'
c get 'GRDPOINT'
c put 'GRDPOINT'
c get 'DIPPOINT'
c put 'DIPPOINT'
c get 'POLPOINT'
c put 'POLPOINT'

#include "flags.h"

      SUBROUTINE UPD_FD(NATOM,DOIT,IMORE,
     &                  ENGPT,GRDPT,DIPPT,POLPT,
     &                  IPTTYPE,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*4 DOIT
      DIMENSION ENGPT(9*NATOM*NATOM),GRDPT(3*NATOM,9*NATOM*NATOM),
     &          DIPPT(3,9*NATOM*NATOM),POLPT(9,9*NATOM*NATOM)
      DIMENSION IPTTYPE(9*NATOM*NATOM),DSCR(NDSCR)

      DIMENSION DIPXYZ(3)
      DIMENSION POLXYZ(3,3)
      LOGICAL PRINTQ

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
#include "control.com"

      irank = 0   ! simulate single-proc execution under ACESII.
      nprocs = 1

      PRINTQ=(IFLAGS(1).GE.10)
      NSIZE=3*NATOM

      CALL IGETREC(20,'JOBARC','NUMPOINT',1,NPOINT)

c   o read point type record and the list of displacements
      CALL IGETREC(20,'JOBARC','FDCALCTP',NPOINT,IPTTYPE)

c Point types can be one of the following:
c  > 0 : this point will be done (1+iRank will be done by this process)
c  = 0 : this point is skipped
c  < 0 : this point is done (-1-iRank was done by this process)

c   o find first entry which must be calculated
      inext = 1
      do while ((ipttype(inext).ne.1+irank).and.(inext.le.npoint))
         inext = inext + 1
      end do

c   o find last entry which was calculated
      ilast = inext-1
      do while ((ipttype(ilast).ne.-1-irank).and.(ilast.gt.0))
         ilast = ilast - 1
      end do

c   o prepare info for the next calculation
      IF (INEXT.NE.NPOINT+1) THEN
         CALL DGETREC(20,'JOBARC','FDCOORDS',NSIZE*INEXT,DSCR)
         ILOC=1+(INEXT-1)*NSIZE
         CALL DPUTREC(20,'JOBARC','NEXTGEOM',NSIZE,DSCR(ILOC))
         IPTTYPE(INEXT)=-IPTTYPE(INEXT)
         CALL IPUTREC(20,'JOBARC','FDCALCTP',INEXT,IPTTYPE)
         IMORE=1
c      o tag the last displacement
         inext = inext+1
         do while ((ipttype(inext).ne.1+irank).and.(inext.le.npoint))
            inext = inext + 1
         end do
         if (inext.eq.npoint+1) then
            CALL IPUTREC(20,'JOBARC','LASTGEOM',1,1)
         end if
      ELSE
         IF (GRADONLY) THEN
            CALL DGETREC(20,'JOBARC','FDCOORDS',NSIZE*ILAST,DSCR)
         END IF
         IMORE=0
      END IF

c   o process info from the last calculation
      IF (ILAST.NE.0) THEN

c      o update energy vector
         CALL DGETREC(20,'JOBARC','ENGPOINT',ILAST,ENGPT)
         IF (IFLAGS(h_IFLAGS_excite).EQ.0) THEN
            CALL DGETREC(20,'JOBARC','TOTENERG',1,ENGPT(ILAST))
         ELSE
            CALL DGETREC(20,'JOBARC','TOTENER2',1,ENGPT(ILAST))
         END IF
         CALL DPUTREC(20,'JOBARC','ENGPOINT',ILAST,ENGPT)

         IF (GRADONLY) THEN

            lFREE = 1

c         o assign the indices, but move the last geom before using DSCR
            lGEOM = lFREE
            lFREE = lFREE + NSIZE
            lGRDXYZ = lFREE
            lFREE = lFREE + NSIZE
            lGRDINT = lFREE
            lFREE = lFREE + NSIZE
            lIMAP = lFREE
            lFREE = lFREE + (NATOM+IINTFP-1)/IINTFP
            NDSCRLFT = NDSCR+1-lFREE

            ILOC=1+(ILAST-1)*NSIZE
            call c_memmove(DSCR(lGEOM),DSCR(ILOC),IFLTLN*NSIZE)

            CALL GETGRD(NATOM,DOIT,DSCR(lGEOM),
     &                  DSCR(lGRDXYZ),DSCR(lGRDINT),DIPXYZ,POLXYZ,
     &                  DSCR(lIMAP),DSCR(lFREE),NDSCRLFT,PRINTQ)

c         o update the FD arrays
            CALL DGETREC(20,'JOBARC','GRDPOINT',NSIZE*ILAST,GRDPT)
            CALL XDCOPY(NSIZE,DSCR(lGRDINT),1,GRDPT(1,ILAST),1)
            CALL DPUTREC(20,'JOBARC','GRDPOINT',NSIZE*ILAST,GRDPT)
            CALL DGETREC(20,'JOBARC','DIPPOINT',3*ILAST,DIPPT)
            CALL XDCOPY(3,DIPXYZ,1,DIPPT(1,ILAST),1)
            CALL DPUTREC(20,'JOBARC','DIPPOINT',3*ILAST,DIPPT)
            CALL DGETREC(20,'JOBARC','POLPOINT',9*ILAST,POLPT)
            CALL XDCOPY(9,POLXYZ,1,POLPT(1,ILAST),1)
            CALL DPUTREC(20,'JOBARC','POLPOINT',9*ILAST,POLPT)

c        END IF (GRADONLY)
         END IF

c     END IF (ILAST.NE.0)
      END IF

      RETURN
      END

