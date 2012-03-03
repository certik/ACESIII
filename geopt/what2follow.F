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

      SUBROUTINE WHAT2FOLLOW(HESMOD, DIAGHES, QSTLST_TANGENT,
     &                       SCRATCH, NX, NOPT, IMODE)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      DOUBLE PRECISION OVRLP

      PARAMETER (OVRLP_MIN = 8.0D-1)
      DATA ZERO /0.0D0/
      DATA TWO  /2.0D0/

      DIMENSION HESMOD(NOPT, NOPT), DIAGHES(NOPT, NOPT),
     &          SCRATCH(NX*NX),  QSTLST_TANGENT(NOPT)
C
C Let's get the eigenvector of the Hessian that has the largest overlap
C with the LST or QST tangent. If there is no eigenvector that has
C an overlap of larger than 10^(-8) or twice as big as the next
C largest overlap, then follow the eigenvector corresponding to
C the smallest eigenvalue. If the user wants to follow something
C other than said eigenvectors, then there is no point of doing
C LST or QST.
C
      ZOVRLP1 = ZERO
      ZOVRLP2 = ZERO

      DO 10 IOPT = 1, NOPT
         OVRLP = DDOT(NOPT,QSTLST_TANGENT,1,DIAGHES(1,IOPT),1)
         IF (DABS(OVRLP).GT.ZOVRLP1) THEN
            ZOVRLP2 = ZOVRLP1
            ZOVRLP1 = DABS(OVRLP)
            MEFOLOW = IOPT
         ELSE IF (DABS(OVRLP) .GT. ZOVRLP2) THEN
            ZOVRLP2 = DABS(OVRLP)
         END IF
 10   CONTINUE
C
c      Write(6,*) "The maximum overlaps = ", ZOVRLP1,ZOVRLP2,
c     &                                      MEFOLOW, OVRLP_MIN
C
      IF (ZOVRLP1 .GT. OVRLP_MIN .OR. ZOVRLP1 .GT. TWO*ZOVRLP2) THEN
         IMODE = MEFOLOW
C
C The pure eigenvector following logic. ie. follow the eiegenvector corresponds
C to the lowest eigenvalue (or user's choice). Find out the eigenvector
C for the lowest eigenvalue (IVEC=1). That is the default.
C
      ELSE

          DO 30 J = 1, NOPT
             SCRATCH(J) = HESMOD(J, J)
 30       CONTINUE
          Z = 100.D0
          DO 50 I = 1, NOPT
             IF (SCRATCH(I) .LT. Z) THEN
                IMODE = I
                Z = SCRATCH(I)
             END IF
 50       CONTINUE

      ENDIF

      RETURN
      END

