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
#include "flags.h"
      SUBROUTINE GEN_APPROX_FREQ(SCRATCH, BMATRIX, TMP, HES_INTACT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "io_units.par"
#include "fnamelen.par"
#include "machsp.com"
#include "mxatms.par"
#include "coord.com"
#include "jodaflags.com"
C
      LOGICAL XYZIN, NWFINDIF
      CHARACTER*1 FLAG(3*MXATMS)
      INTEGER TOTREDNCO
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
      DOUBLE PRECISION BMATRIX(NXM6*3*NATOMS), TMP(NXM6*3*NATOMS),
     &                 SCRATCH(9*NATOMS*NATOMS),
     &                 HES_INTACT(3*NATOMS, 3*NATOMS)

C Transform the internal coordinate Hessian to the Cartesian Hessian in
C order to carry out mass weighing (it is nearly impossible to mass
C weigh internal Hessians). Notice the assumption made in the
C transformation: the gradient is zero, hence no derivative of the B
C matrix is required. Since this is done only at the stationary points,
C the assumption (or the condition) is satisfied and the transformation
C is valid. Ajith Perera 05/2005.

CSSS      CALL OUTPUT(HES_INTACT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
      IF (Iflags2(h_IFLAGS2_geom_opt) .GE. 3) THEN
C
C It is more than safer to always use 3*NATOMS explicitly when
C redundent internal optimizations are performed. Ajith perera,
C 01/2007.

         CALL IGETREC(20, 'JOBARC','REDNCORD', 1, TOTREDNCO)
         LENGTH_BMAT = 3*NATOMS*TOTREDNCO
         CALL DGETREC(-1, 'JOBARC', 'BMATRIXT', LENGTH_BMAT,
     &              TMP)
         CALL ZERO(SCRATCH, 9*NATOMS*NATOMS)
         CALL XGEMM("N", "N", 3*NATOMS, TOTREDNCO, TOTREDNCO, 1.0D0,
     &              TMP, 3*NATOMS, HES_INTACT, TOTREDNCO, 0.0D0,
     &              SCRATCH, 3*NATOMS)
         CALL TRANSP(TMP, BMATRIX, TOTREDNCO, 3*NATOMS)
         CALL ZERO(HES_INTACT, 9*NATOMS*NATOMS)
         CALL XGEMM("N", "N", 3*NATOMS, 3*NATOMS, TOTREDNCO, 1.0D0, 
     &               SCRATCH, 3*NATOMS, BMATRIX, TOTREDNCO, 0.0D0, 
     &               HES_INTACT, 3*NATOMS)
C
      ELSE
C
         CALL DGETREC(20, 'JOBARC', 'BMATRIX ', NXM6*NX,
     &               BMATRIX)
CSSS         CALL OUTPUT(BMATRIX, 1, NXM6, 1, NX, NXM6, NX, 1)
C
         CALL TRANSP(BMATRIX, TMP, NX, NXM6)
         CALL ZERO(SCRATCH, NX*NX)
         CALL XGEMM("N", "N", NX, NXM6, NXM6, 1.0D0, TMP, NX,
     &               HES_INTACT, NXM6, 0.0D0, SCRATCH, NX)
         CALL ZERO(HES_INTACT, NX*NX)
         CALL XGEMM("N", "N", NX, NX, NXM6, 1.0D0, SCRATCH, NX,
     &               BMATRIX, NXM6, 0.0D0, HES_INTACT, NX)
C
CSSS         CALL OUTPUT(HES_INTACT, 1, NX, 1, NX, NX, NX, 1)
      ENDIF
C
C Mass weight the Cartesian Hessian
C
      DO I = 1, 3*NATOMS
         DO J = 1, 3*NATOMS
            FACT = DSQRT(ATMASS(1+(I-1)/3)*ATMASS(1+(J-1)/3))
            IF (FACT .LT. 1.0D-3)THEN
                HES_INTACT(I,J) = 0.D0
            ELSE
                HES_INTACT(I,J) = HES_INTACT(I,J)/FACT
            ENDIF
         ENDDO
      ENDDO
C
C Diagonalize the mass-weighted Hessian and tag the negative
C eigenvalues with "i" to indicate an imaginary frequency and
C convert them to the correspnding frequencies in cm-1 units.
C
      CALL EIG(HES_INTACT, SCRATCH, 3*NATOMS, 3*NATOMS, 0)
C
      DO I = 1, 3*NATOMS
         FLAG(I) = " "
         IF (HES_INTACT(I,I) .LT. 0.0D0) FLAG(I) = "i"
         HES_INTACT(I,I) = DSQRT(DABS(HES_INTACT(I,I)))*5.14048D03
         SCRATCH(I) = HES_INTACT(I,I)
      END DO
C
C Print the frequencies along with a warning message.
C
      WRITE(LUOUT,10)
   10 FORMAT(3X, "Frequencies of the updated Hessian at convergence:")
      WRITE(LUOUT,*)
C
      WRITE(LUOUT,'(11(F6.0,A))')
     &     (SCRATCH(IMODE), FLAG(IMODE), IMODE=1,3*NATOMS)
      WRITE(LUOUT,*)
      WRITE(LUOUT,20)
   20 FORMAT(3X, "Warning: Frequencies based on the updated Hessian",
     &" are not very accurate.",/,12X, "A separate frequency",
     &" calculation is required to get correct values.")
      WRITE(LUOUT,*)
C
      RETURN
      END

