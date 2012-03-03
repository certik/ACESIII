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
c get ROTREC(IXYZ)
c put ROTREC(IXYZ)
c get 'COMPNSYQ'
c get 'COMPSYMQ'
c put 'COMPSYMQ'

      SUBROUTINE ROTPRJ(NATOM,TPROJ,SCR,SYMQ,ATMASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SCR(*),SYMQ(*),TPROJ(*),ATMASS(*)

      CHARACTER*8 ROTREC(3)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA ONE /1.0/
      DATA TOL /1.D-8/
      DATA ROTREC /'ROTVECX ','ROTVECY ','ROTVECZ '/

      CALL IGETREC(20,'JOBARC','LINEAR  ',1,ILINEAR)

C CONSTRUCT THE ROTATIONAL PROJECTOR = 1 - SUM |Rq><Rq|
C                                          xyz

      NSIZE=3*NATOM
      CALL ZERO(TPROJ,NSIZE*NSIZE)
      CALL DGETREC(20,'JOBARC','ATOMMASS',NATOM,ATMASS)
      IOFF=1
      DO I=1,NSIZE
         TPROJ(IOFF)=ONE
         IOFF=IOFF+NSIZE+1
      END DO
      DO IXYZ=1,3-ILINEAR
c      o mass-weight and write back to disk
         CALL DGETREC(20,'JOBARC',ROTREC(IXYZ),NSIZE,SCR)
         IOFF=1
         DO I=1,NATOM
            Z=SQRT(ATMASS(I))
            CALL XDSCAL(3,Z,SCR(IOFF),1)
            IOFF=IOFF+3
         END DO
         Z=DNRM2(NSIZE,SCR,1)
         FACT=ONE/Z
         CALL XDSCAL(NSIZE,FACT,SCR,1)
         CALL DPUTREC(20,'JOBARC',ROTREC(IXYZ),NSIZE,SCR)
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
         END IF
         IOFF=IOFF+NSIZE
      END DO

      CALL DPUTREC(20,'JOBARC','COMPSYMQ',NSIZE*NSYMOLD,SYMQ)

      RETURN
      END

