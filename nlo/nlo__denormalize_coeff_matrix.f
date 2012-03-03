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
         SUBROUTINE  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +                    ( NBAS,
     +                      DENORM,
     +
     +                              COEFFS )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__DENORMALIZE_COEFF_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine denormalizes the coefficient matrix such
C                that the natural orbitals are expressed in the original
C                AO basis set. This is necessary sometimes when the
C                overlap and density matrices in the original AO basis
C                set have both been normalized and the natural orbital
C                coefficients are based on the normalized AO's.
C
C                The denormalized coefficient matrix is simply obtained
C                by multiplying each row with the appropriate
C                denormalization factor (which is equal to the square
C                root of the diagonal unnormalized overlap matrix).
C
C
C                Input:
C
C                     NBAS = size of AO basis.
C                   DENORM = denormalization coefficients.
C                   COEFFS = original coefficient matrix.
C
C                Output:
C
C                   COEFFS = denormalized coefficient matrix.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    I,J
         INTEGER    NBAS

         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  X

         DOUBLE PRECISION  DENORM (1:NBAS)

         DOUBLE PRECISION  COEFFS (1:NBAS,1:NBAS)

         PARAMETER    (ONE  =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform denormalization. Skip multiplication of
C                rows if denormalization factor is equal to 1.
C
C
         DO 10 I = 1,NBAS
            X = DENORM (I)
            IF (X.NE.ONE) THEN
                DO 20 J = 1,NBAS
                   COEFFS (I,J) = (ONE/X) * COEFFS (I,J)
   20           CONTINUE
            END IF
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
