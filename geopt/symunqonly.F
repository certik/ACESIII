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
C A fairly simple routine to identify symmetry unique values in the
C array RIN and put the symmetry unique elements to the begining
C of the RIN array followed by locating the elements in the RIN
C array that are being optimized and putting them in ROUT array
C in the proper (to ACES II) order.  Ajith Perera, 01/2005.
C
      SUBROUTINE SYMUNQONLY(RIN, ROUT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "coord.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT

      DIMENSION RIN(NX), ROUT(NX)

      DO I = 1, NUNIQUE
         IP = IUNIQUE(I)
         Z  = RIN(IP)
         DO J = 1, NEQ(IP)
            Z = Z + RIN(IEQUIV(I,J))
         END DO
         FIAVE  = Z/(NEQ(IP)+1)
         RIN(I) = FIAVE
      END DO
      DO IUNQ = 1, NUNIQUE
         ROUT(IUNQ) = RIN(NOPTI(IUNQ))
      END DO
      RETURN
      END

