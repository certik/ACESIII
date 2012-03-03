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

      SUBROUTINE GET_MOPAC_HESS(A, SCRATCH, IMAP, HC, HI, AT, NATOMS,
     &                          NX, NXM6, NREAL, NREAL3, NREAL3M6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "machsp.com"
#include "fnamelen.par"
#include "io_units.par"

      DOUBLE PRECISION HC(NX*NX),SCRATCH(NX*NX), A(NX*NXM6),
     &                 AT(NXM6*NX), HI(NXM6, NXM6)
      INTEGER IMAP(NATOMS)
      COMMON /OPTCTL/ IPRNT, INR, IVEC, IDIE, ICURVY, IMXSTP, ISTCRT,
     &                IVIB, ICONTL, IRECAL, INTTYP, IDISFD, IGRDFD,
     &                ICNTYP, ISYM, IBASIS, XYZTol
      CHARACTER*(fnamelen) FNAME
      PARAMETER (IMOPAC_UNIT = 99)

      CALL GFNAME('NDDO_HES',FNAME, ILENGTH)
      OPEN(UNIT=IMOPAC_UNIT, FILE=FNAME(1:ILENGTH), FORM=
     &     "UNFORMATTED", STATUS="OLD")

      LTRIN = NREAL*3 *(NREAL*3 + 1)/2
      READ(IMOPAC_UNIT) (SCRATCH(I), I=1, LTRIN)
      CALL EXPND2(SCRATCH, HC, NREAL3)
c      Print*, "The Hessian just read from MOPAC"
CSSS      CALL OUTPUT(SCRATCH, 1, NREAL3, 1, NREAL3, NREAL3,
CSSS     &            NREAL3, 1)

C Take care of the cases that invlove dummy or ghost atoms.
      NDUMMYN = NX*NX
      IF (NREAL3.NE.NX) THEN
         NREALN  = NREAL3*NREAL3
         CALL ZERO(SCRATCH, NDUMMYN)
         CALL IGETREC(1,'JOBARC','MAP2ZMAT',NATOMS,IMAP)
         IOFF = 1
         DO IATREL = 1, NATOMS
            IATMZMAT = IMAP(IATREL)
            DO IXYZ = 1, 3
               ICOL = IXYZ + (IATMZMAT-1)*3
               DO JATREL = 1, NATOMS
                  JATMZMAT = IMAP(JATREL)
                  IF (IATMZMAT .NE.0 .AND. JATMZMAT .NE. 0) THEN
                     IROW = 1 + (JATMZMAT-1)*3
                     CALL BLKCPY(HC(IOFF), 3, 1, SCRATCH, NX,
     &                           NX, IROW, ICOL)
                     IOFF = IOFF + 3
                  END IF
               END DO
            END DO
         END DO
         CALL DCOPY(NDUMMYN, SCRATCH, 1, HC, 1)
      END IF

C According to Keith Runge, the MOPAC Hessian is in mdynes/picometer.
C The conversion factor was provided by him.  Priv. comm.  12/2004
C
      CONVERT_TO_AU = (6.2414D0/27.2113957D0)*
     &                (1.0D0/0.529177249D0)**2
      CALL DSCAL(NDUMMYN, CONVERT_TO_AU, HC, 1)
C
C Convert to internal coordinates using the A matrix.
C
      CALL MODMATMUL(SCRATCH,HC,A,NX,NX,NXM6,NX,NX,NXM6)
      CALL MTRANSP(A,AT,NX,NXM6,NX,NXM6)
      CALL MODMATMUL(HI,AT,SCRATCH,NXM6,NX,NXM6,NXM6,NX,NXM6)

      IF (IPRNT.GE.10) THEN
         WRITE(LUOUT,110)
 110     FORMAT(T3,' @GET_MOPAC_HESS, Full Cartesian Hessian: ')
         CALL OUTPUT(HC, 1, NX, 1, NX, NX, NX, 1)
         Write(Luout, 111)
 111     FORMAT(T3,' @GET_MOPAC_HESS, Full internal coordinate
     & Hessian: ')
          CALL OUTPUT(HI, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
      END IF
      CLOSE(IMOPAC_UNIT)

      RETURN
      END

