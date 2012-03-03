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
c get 'COMPNSYQ'
c get 'COMPSYMQ'
c put 'COMPSYMQ'

      SUBROUTINE TRAPRJ(NATOM,TPROJ,SCR,SYMQ,ATMASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SCR(1),SYMQ(1),TPROJ(1),ATMASS(1)

      CHARACTER*8 TRAREC(3)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA ONE /1.0/
      DATA TOL /1.D-10/
      DATA TRAREC /'TRAVECX ','TRAVECY ','TRAVECZ '/

      CALL DGETREC(20,'JOBARC','ATOMMASS',NATOM,ATMASS)

C CONSTRUCT THE TRANSLATIONAL PROJECTOR IN MASS WEIGHTED COORDS
C           1 - SUM |T><T|

      NSIZE=3*NATOM
      CALL ZERO(TPROJ,NSIZE*NSIZE)
      IOFF=1
      DO I=1,NSIZE
         TPROJ(IOFF)=ONE
         IOFF=IOFF+NSIZE+1
      END DO
      DO IXYZ=1,3
         CALL ZERO(SCR,NSIZE)
         ndx = IXYZ
         DO I=1,NATOM
            SCR(ndx)=SQRT(ATMASS(I))
            ndx = ndx + 3
         END DO
         FACT=ONE/DNRM2(NSIZE,SCR,1)
         CALL XDSCAL(NSIZE,FACT,SCR,1)
         CALL DPUTREC(20,'JOBARC',TRAREC(IXYZ),NSIZE,SCR)
         CALL XGEMM('N','N',NSIZE,NSIZE,1,
     &              -1.d0,SCR,  NSIZE,
     &                    SCR,  1,
     &              ONE,  TPROJ,NSIZE)
      END DO

c   o hit symmetry adapted coordinates with projector
      CALL IGETREC(20,'JOBARC','COMPNSYQ',1,NSYMOLD)
      CALL DGETREC(20,'JOBARC','COMPSYMQ',NSYMOLD*NSIZE,SCR)
      CALL XGEMM('N','N',NSIZE,NSYMOLD,NSIZE,ONE,TPROJ,NSIZE,
     &           SCR,NSIZE,0.d0,SYMQ,NSIZE)

c   o renormalize
      IOFF=1
      DO I=1,NSYMOLD
         X=DNRM2(NSIZE,SYMQ(IOFF),1)
         IF (ABS(X).GT.TOL) THEN
            Z=ONE/X
            CALL XDSCAL(NSIZE,Z,SYMQ(IOFF),1)
            IOFF=IOFF+NSIZE
         END IF
      END DO

      CALL DPUTREC(20,'JOBARC','COMPSYMQ',NSYMOLD*NSIZE,SYMQ)

      RETURN
      END

