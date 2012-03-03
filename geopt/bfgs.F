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

c The Broyden-Fletcher-Goldfarb-Shanno update.

      SUBROUTINE BFGS(V,H,SCRATCH,STEP,TBT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "io_units.par"
#include "fnamelen.par"
#include "jodaflags.com"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION V(NXM6),H(NXM6*NXM6),SCRATCH(NX*NX),TBT(3*NXM6*NXM6),
     $          STEP(NXM6)

      N2=1+NXM6
      N3=1+2*NXM6
      N4=1+3*NXM6
      NX2=1+NXM6*NXM6
C
C SCRATCH(N2) = DG
C
      CALL VADD(SCRATCH(N2),V,SCRATCH(1),NXM6,-1.D0)
C
C TBT(1) =  [1/DG{dot}DX][DGxDG]
C
      CALL MODMATMUL(TBT(1),SCRATCH(N2),SCRATCH(N2),NXM6,1,NXM6,
     $               NXM6,1,NXM6)
      X0=1.D0/xdot(NXM6,SCRATCH(N2),1,STEP,1)
      Write(6,*) "IN BFGS X0", X0
      CALL xscal(NXM6*NXM6,X0,TBT(1),1)
C
C SCRATCH(N3) = HDX, SCRATCH(N4) = DXH, TBT(NX2)= HDX{x}DXH
C
      CALL MODMATMUL(SCRATCH(N3),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
      CALL MODMATMUL(SCRATCH(N4),STEP,H,1,NXM6,NXM6,1,NXM6,NXM6)
      CALL MODMATMUL(TBT(NX2),SCRATCH(N3),SCRATCH(N4),NXM6,1,NXM6,
     $               NXM6,1,NXM6)
C
C Z0 = 1/DXHDX, TBT(NX2) = [1/DXHDX][HDX{x}DXH]
C
      Z0=1.D0/xdot(NXM6,STEP,1,SCRATCH(N3),1)
      Write(6,*) "IN BFGS Z0", Z0
      CALL xscal(NXM6*NXM6,Z0,TBT(NX2),1)
C
C H(updated) = H(old) + [1/DG{dot}DX][DGxDG] - [1/DXHDX][HDX{x}DXH]
C
      CALL VADD(TBT(1),TBT(1),TBT(NX2),NXM6*NXM6,-1.D0)
      CALL VADD(H,H,TBT(1),NXM6*NXM6,1.D0)
C
C DON'T WRITE ENTIRE HESSIAN UNLESS SPECIFICALLY REQUESTED.
C FOR NOW, JUST USE AN IN-CODE PARAMETER.
C
#ifdef _DEBUG_LVLM1
         Print*, "-----The BFGS updated Hessian-----"
         CALL HESSOUT(H,NXM6,NXM6,0)
#endif
      RETURN
      END

