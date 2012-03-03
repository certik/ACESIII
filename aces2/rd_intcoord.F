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

#include "iachar.h"
#include "flags.h"

      subroutine rd_intcoord(rint)
      implicit none

      DOUBLE PRECISION DACOS
      INTRINSIC DACOS
      INTEGER LINBLNK
      EXTERNAL LINBLNK

#include "io_units.par" /* for LUZ */
#include "mxatms.par"
#include "linelen.par"
#include "cbchar.com"
#include "coord.com"
C
      DOUBLE PRECISION dSmallNumber
      DOUBLE PRECISION RINT(3*MXATMS)
      PARAMETER (dSmallNumber=5.361d-21)

      INTEGER        NX,NXM6,IARCH,NCYCLE,NUNIQUE,NOPT
      COMMON /USINT/ NX,NXM6,IARCH,NCYCLE,NUNIQUE,NOPT

      INTEGER        IFLAGS(100),IFLAGS2(500)
      COMMON /FLAGS/ IFLAGS,IFLAGS2

      DOUBLE PRECISION RTOA, dTmp
      INTEGER I, J, IBOMB, ISTP
      CHARACTER*(linelen) DUMSTR
      CHARACTER*5 CHTEST
      CHARACTER*1 czTab
      LOGICAL PRINTI

c ----------------------------------------------------------------------

C NEW (5/89) VARIABLE INPUT. NO LONGER HAVE TO INCLUDE REDUNDANCIES, AN
C UNIQUE COORDINATES CAN BE READ IN ANY ORDER. VARIABLE NAMES AND VAL
C ARE SEPARATED BY AN EQUALS (=) SIGN. MUCH BETTER THAN THE WAY IT WA
C DONE BEFORE. START BY FILLING R VECTOR WITH BOGUS VALUES.

      IF (IFLAGS(h_IFLAGS_coordinates).EQ.2) THEN

c Ajith Perera 12/2001.
C This block corresponds to the COORINATES=XYZ2INT option. To
C understand what this means, one needs to know the difference
C between COORDINATES=CARTESIAN and COORDINATES=XYZINT. In the first
C case, our input is the atomic symbol followed by the Cartesian
C coordinates of the each atom. There is no information about
C the internal coordinates. In the latter case, as one would
C have done for entering internal coordinates, one would enter
C the connectivity data in ZMAT format followed by (one blank
C in between) the Cartesian coordinates of the each atom in free
C format form.
C At present, this option does not understand the symmetry unless
C the user imposes it by ZMAT input. The advantage of having
C this capability is that one can use the Cartesian coordinates
C directly without first converting them to internal.
C In the future, we plan to extend this to include symmetry
C without user specfications. That way one can do optimiations
C with Cartesian input.
C
C Originaly, we had a routine called XTOR2.F which was pretty
C much doing what XTOR.F was doing except for the conversion of
C angles to radians. I obsoleted the XTOR2.F and extended the
C XTOR.F. The IPRT flag (don't know why it was named that)
C controls whether the conversion is needed (IPRT==1) or not
C (IPRT!=1).
         WRITE(LUOUT,*) '@RD_INTCOORD: Input coordinates are ',
     &      'transformed to internal representation.'
         READ(LUZ,*) (Q(I),I=1,NX)
         CALL XTOR(R,1)

      ELSE

         CALL IGETREC(1,'JOBARC','FIRSTRUN',1,I)
         PRINTI = (I.NE.0)
         czTab  = achar(_IACHAR_TAB)
         RTOA   = 180.D0/DACOS(-1.D0)
         IBOMB  = 0
         DO I = 1, NX
            R(I) = dSmallNumber
         END DO

         READ(LUZ,'(A)',END=987) DUMSTR
         do while (linblnk(dumstr(1:80)).ne.0)
            ISTP = INDEX(DUMSTR(1:),'=')
            I    = INDEX(DUMSTR(1:),'*')
            if ((i.ne.0).or.(istp.eq.0)) then
               write(*,*) '@RD_INTCOORD: Coordinate list is not ',
     &                    'terminated with a blank line.'
               call errex
            end if
            i = 1
            do while ((i.lt.istp).and.
     &                (dumstr(i:i).eq.' '.or.dumstr(i:i).eq.czTab))
               i = i + 1
            end do
            if (i.eq.istp) then
               write(*,*) '@RD_INTCOORD: Coordinate name is missing.'
               call errex
            end if
            j = i + 1
            do while ((j.lt.istp).and.
     &                (dumstr(j:j).ne.' '.and.dumstr(j:j).ne.czTab))
               j = j + 1
            end do
            CHTEST = DUMSTR(I:J-1)
            READ(DUMSTR(ISTP+1:),*) dTmp
            ISTP = 0
            I = 1
            do while (ISTP.eq.0.and.I.le.NXM6)
               IF (CHTEST.EQ.VARNAM(ISQUASH(I))) THEN
                  ISTP = 1
                  IF (RINT(I).NE.dSmallNumber) THEN
                     IBOMB = 1
                     WRITE(*,*) '@RD_INTCOORD: ',CHTEST,
     &                          ' multiply defined.'
                  ELSE IF (VARNAM(ISQUASH(I)).EQ.'TDA  ') THEN
                     RINT(I) = RTOA*DACOS(-1.0D0/3.0D0)
                     IF (PRINTI) WRITE(LUOUT,*)
     & '@RD_INTCOORD: Variable TDA set to exact tetrahedral angle.'
                  ELSE IF (VARNAM(ISQUASH(I)).EQ.'IHA  ') THEN
                     RINT(I) = RTOA*DACOS(1.D0/DSQRT(5.D0))
                     IF (PRINTI) WRITE(LUOUT,*)
     & '@RD_INTCOORD: Variable IHA set to exact icosahedral angle.'
                  ELSE IF (VARNAM(ISQUASH(I)).EQ.'OA   ') THEN
                     RINT(I) = RTOA*DACOS(-1.D0/3.D0)*0.5D0
                     IF (PRINTI) WRITE(LUOUT,*)
     & '@RD_INTCOORD: Variable OA set to 0.5 * exact tetrahedral angle.'
                  ELSE
                     R(I) = dTmp
                  END IF
               END IF
               I = I + 1
c           end do while (ISTP.eq.0.and.I.le.NXM6)
            end do
            IF (ISTP.EQ.0) THEN
               IBOMB = 1
               WRITE(*,*) '@RD_INTCOORD: ',CHTEST,' not in Z-matrix.'
            END IF
            READ(LUZ,'(A)',END=987) DUMSTR
c        end do while (linblnk(dumstr(1:80)).ne.0)
         end do
  987    CONTINUE

c      o copy unique coordinate definitions
         DO I = 1, NUNIQUE
            DO J = 1, NEQ(IUNIQUE(I))
               RINT(IEQUIV(I,J)) = RINT(IUNIQUE(I))
            END DO
         END DO

c      o check if all coordinates have been defined
         DO I = 1, NXM6
            IF (RINT(I).EQ.dSmallNumber) THEN
               IBOMB = 1
               WRITE(*,*)
     &               '@RD_INTCOORD: ',VARNAM(ISQUASH(I)),' is not ',
     &               'defined at ',ISQUASH(I)
            END IF
         END DO

c      o stop on all errors
         IF (IBOMB.EQ.1) CALL ERREX

      END IF

c   o unsquash R
      CALL USQUSH(R,NXM6)

      return
      end

