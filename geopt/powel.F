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

      SUBROUTINE POWEL(V,H,SCRATCH,STEP,TBT, PHI)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "io_units.par"
#include "fnamelen.par"
#include "jodaflags.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT

      DIMENSION V(NXM6),H(NXM6*NXM6),SCRATCH(NX*NX),TBT(3*NXM6*NXM6),
     $   STEP(NXM6)
C
C X = 1/DX{dot}DX where DX is the step size.
C
      X=1.D0/xdot(NXM6,STEP,1,STEP,1)

      N2=NXM6+1
      N3=2*NXM6+1
      N4=3*NXM6+1
C
C SCRATCH(N2) = HDX where H is the Hessian, SCRATCH(N3) = DG (gradient
C difference), SCRATCH(N4) = (DG - HDX)
C
      CALL MODMATMUL(SCRATCH(N2),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
      CALL VADD(SCRATCH(N3),V,SCRATCH(1),NXM6,-1.D0)
      CALL VADD(SCRATCH(N4),SCRATCH(N3),SCRATCH(N2),NXM6,-1.D0)

      N2=NXM6*NXM6+1
      N3=NXM6*NXM6+N2
C
C TBT(1) = DXDX, TBT(N2) = (DG - HDX)DX
C
      CALL MODMATMUL(TBT(1),STEP,STEP,NXM6,1,NXM6,NXM6,1,NXM6)
      CALL MODMATMUL(TBT(N2),SCRATCH(N4),STEP,NXM6,1,NXM6,NXM6,1,
     &               NXM6)
C
C Take the transpose of (DG - HDX)DX, TBT(N3) = [(DG - HDX)DX]^{t}
C
      CALL MTRANSP(TBT(N2),TBT(N3),NXM6,NXM6,NXM6,NXM6)
C
C Built the (DG - HDX)dotDX and scale the DXDX, TBT(1) =
c [DXdot(DG - HDX)]DXDX
C
      dtmp = xdot(NXM6,SCRATCH(N4),1,STEP,1)
      CALL xscal(NXM6*NXM6,dtmp*X,TBT(1),1)
C
C TBT(N3) = (DG - HDX)DX + [(DG - HDX)DX]^{t} - [DXdot(DG - HDX)]DXDX/DX{dot}DX
C
      CALL VADD(TBT(N3),TBT(N3),TBT(N2),NXM6*NXM6,1.D0)
      CALL VADD(TBT(N3),TBT(N3),TBT(1),NXM6*NXM6,-1.D0)
C
C H(updated) = H(old) + [(DG - HDX)DX + [(DG - HDX)DX]^{t}]/DX{dot}DX +
C              DXdot(DG - HDX)]DXDX/[DX{dot}DX]^{2}
C
      CALL xscal(NXM6*NXM6,X,TBT(N3),1)
      CALL XDAXPY(NXM6*NXM6, PHI, TBT(N3), 1, H, 1)
C
C DON'T WRITE ENTIRE HESSIAN UNLESS SPECIFICALLY REQUESTED.  FOR
C NOW, JUST USE AN IN-CODE PARAMETER.
C
      IF (IFLAGS(1).GE.10) THEN
         Print*, "-----The Powel updated Hessian-----"
         CALL HESSOUT(H,NXM6,NXM6,0)
      END IF

      RETURN
      END

