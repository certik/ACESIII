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

C THIS ROUTINE TESTS ALL SYMMETRY OPERATIONS OF THE GROUP ON A
C PARTICULAR SET OF COORDINATES. ALSO RETURNS POINTER LISTS
C FOR ALL SYMMETRY OPERATIONS.

#define _DEBUG_TSTOPS

      SUBROUTINE TSTOPS(NORDER,Q,V,CLSTYP,IPTR,SCR,NORD,NATOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates

      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /TOLERS/ SYMTOL,DEGTOL
 
#include "io_units.par"
 
      DIMENSION V(9*NORDER),Q(3*NATOM),SCR(3*NATOM)
      DIMENSION NORD(2*NATOM),IPTR(NATOM,NORDER)
      INTEGER CLSTYP(NORDER)
      CHARACTER*6 RESULT

      iprnt = 101
      ICOMP = 0
      CALL ZERO(SCR,3*NATOM)
c      CALL IZERO(NORD,NATOM)
    
      DO I = 1, NATOM
         NORD(I) = I
      END DO
      CALL SORTXYZ(Q,SCR,NORD,NATOM)

      IF (IPRNT.GE.100) THEN
         WRITE(LUOUT,999)
         WRITE(LUOUT,2000)
         WRITE(LUOUT,1000)
         WRITE(LUOUT,2000)
      END IF
 999  FORMAT(T3,'@TSTOPS: Results of symmetry operation testing:')
 2000 FORMAT(80('-'))
 1000 FORMAT(T3,' Sym. Op.',T20,' Class ',T30,' Trace ',T40,' Result ')

      ISUCC = 0
      NFAIL = 0
      DO I = 1, NORDER
         itmp = -8 + ( 9 * I )
         CALL XGEMM('T','N',3,NATOM,3,
     &              1.0d0,V(itmp),       3,
     &                    Q,             3,
     &              0.0d0,SCR(3*NATOM+1),3)
         CALL COMPARE2(SCR(3*NATOM+1),SCR,NORD,ICOMP,SYMTOL)
         CALL STPTR(NATOM,NORD,NORD(NATOM+1),IPTR(1,I))
         IF (ICOMP.EQ.0) THEN
            RESULT = 'Passed'
            ISUCC  = ISUCC + 1
         ELSE
            RESULT = 'Failed'
            NFAIL  = NFAIL + 1
#ifdef _DEBUG_TSTOPS
            write(*,*) '@TSTOPS: ',icomp,' coords are above ',symtol
#endif /* _DEBUG_TSTOPS */
         END IF
         IF (IPRNT.GE.100) THEN
            itmp = -8 + ( 9 * I )
            TRACE = V(itmp) + V(itmp+4) + V(itmp+8)
            WRITE(LUOUT,1001) I,CLSTYP(I),TRACE,RESULT
         END IF
 1001    FORMAT(T6,I3,T22,I3,T30,F8.5,T41,A6)
      END DO
      IF (NFAIL.GT.0) THEN
         WRITE(*,*) '@TSTOPS: ',NFAIL,' operations failed.'
         CALL ERREX
      END IF
      IF (IPRNT.GE.100) WRITE(LUOUT,2000)
      IF (IPRNT.GT.5) THEN
         WRITE(*,*) '@TSTOPS: ',ISUCC,' operations verified.'
      END IF

      RETURN
      END

