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
C  Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
C 
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C  USA
      SUBROUTINE SETUP(HES, GRD, HESMOD, GRDMOD, IQFIX, STPMAX, LUOUT)
C
C Build the Hessian and Gradient vectors symmetry coordintes (
C only those internal that are optimized is taking in to account)
C Also do all other initializations.
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "coord.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT, INR, IVEC, IDIE, ICURVY, IMXSTP, ISTCRT, 
     &                IVIB, ICONTL, IRECAL, INTTYP, IDISFD, IGRDFD, 
     &                ICNTYP, ISYM, IBASIS, XYZTol

C
      DIMENSION HES(NXM6, NXM6), GRD(NXM6), HESMOD(NOPT, NOPT),
     &          GRDMOD(NOPT), IQFIX(3*NATOMS, 3*NATOMS)
C
C Let's update the number of cycles right away!
C
      NCYCLE = NCYCLE + 1 
C
      DO 311 I=1, NUNIQUE
         DO 312 J=1, NEQ(IUNIQUE(I))
            IQFIX(IUNIQUE(I), J)=IEQUIV(I,J)
 312     CONTINUE
 311  CONTINUE
C  
      DO 19 I=1, NOPT
         DO 20 J=1, NOPT
C
            HESMOD(I,J)=HES(NOPTI(I), NOPTI(J))
C
            DO 29 K=1, NEQ(NOPTI(I))
               HESMOD(I,J)=HESMOD(I,J)+HES(IQFIX(NOPTI(I),K),NOPTI(J))
 29         CONTINUE
C
            DO 39 L=1, NEQ(NOPTI(J))
               HESMOD(I,J)=HESMOD(I,J)+HES(NOPTI(I),IQFIX(NOPTI(J),L))
 39         CONTINUE
C
            DO 49 K=1, NEQ(NOPTI(I))
               DO 50 L=1, NEQ(NOPTI(J))
                  HESMOD(I, J)=HESMOD(I, J) + HES(IQFIX(NOPTI(I), K),
     &                        IQFIX(NOPTI(J), L))
 50            CONTINUE
 49         CONTINUE
C
            ZD=DSQRT(DFLOAT((NEQ(NOPTI(I))+1)*(NEQ(NOPTI(J))+1)))
            HESMOD(I,J)=HESMOD(I,J)/ZD
C
 20      CONTINUE
 19   CONTINUE
C
      IF (IPRNT .GE. 10) THEN
          Write(6,*)
          WRITE(LUOUT,115)
 115      FORMAT(T3,' Hessian matrix in totally symmetric symmetry ',
     &      'coordinates: ')
          IF(ICURVY.EQ.1) WRITE(LuOut,117)
 117      FORMAT(T3,'            (Curvilinear representation) ')
          CALL HESSOUT(HESMOD, NOPT, NOPT, 1)
      ENDIF
C     
C Form modified Gradient vector and symmetrize it. 
C
      DO 1177 I=1, NOPT
         GRDMOD(I)=DSQRT(DFLOAT(NEQ(NOPTI(I))+1))*GRD(NOPTI(I))
 1177 CONTINUE
C
C The maximum step size can be controlled by the user by setting
C the flag MAX_STEP. The default value for maximum is set 300.
C IMXSTP is first set in mkvmol.F (strange!). We need to move
C these things into a one place that is visible to the developers.
      IF (IMXSTP .EQ. 0) IMXSTP = 300
      STPMAX = DBLE(IMXSTP)/1000.D0
C
      RETURN
      END




