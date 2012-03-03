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

C  *******************************************************
C  Conversion of the hessian from cartesian to red. internal
C  coordinates.
C  Luis Galiano, 07/07/03
C  *******************************************************

      SUBROUTINE CART2INT_HESS(FI,HC,TOTREDNCO,NRATMS,HI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER TOTREDNCO
#include "io_units.par"
#include "fnamelen.par"
#include "mxatms.par"

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "ric_heap.com"

C The following arrays are need to be dynamically allocated
C and deallocated.
C
C DIFTEMP: The is used in COMPKMAT as a temporary array to keep
C          an intermediate. It is of length 9*MXATMS*MXATMS.
C TEMP1  : Antoher temporary array to keep intermediates during
C          the transformations of the Cartesian Hessian. It is
C          of length 3*MAXREDUNCO*MXATMS.
C GMINBT : Keep the G-mat for Hessian transformations. It is of
C          length 9*MXATMS*MXATMS.
C GMINB  : To keep the  transpose of the G-mat, and it is of
C          the same lengthe as GMINBT.
C PFI    : Keep the gradients during the gradient projection.
C          It is of length MAXREDUNCO.
      DIMENSION DIFTEMP(9*NRATMS*NRATMS),FI(TOTREDNCO),PFI(TOTREDNCO),
     &          HC(3*NRATMS,3*NRATMS),
     &          HI(TOTREDNCO*TOTREDNCO),
     &          GMINBT(3*NRATMS*MAXREDUNCO),
     &          TEMP1(TOTREDNCO*3*NRATMS),
     &          GMINB(TOTREDNCO*TOTREDNCO)

      CALL IGETREC(-1,'JOBARC','REDNCORD',1,TOTREDNCO)
      IPRNT=600

C Build the K-matrix and the (H-K) matrix as shown in JCP, 117, 9160, 2002.
#ifdef _DEBUG_LVLM1
      Print*, "The Cartesian force constant matrix"
      CALL OUTPUT(HC,1,3*NRATMS,1,3*NRATMS,3*NRATMS,3*NRATMS, 1)
#endif
      CALL COMPKMAT(dRICHeap(z_DerBMat),NRATMS,TOTREDNCO,FI,HC,DIFTEMP)

      LENGTH_BGMAT =3*NRATMS*TOTREDNCO
      CALL DGETREC(20,'JOBARC','BTGMIN',LENGTH_BGMAT,
     &            dRICHeap(z_BTGInv))
#ifdef _DEBUG_LVL0
         print *, 'The BTGINV matrix'
         CALL OUTPUT(dRICHeap(z_BTGInv), 1, 3*NRATMS, 1, TOTREDNCO, 
     &               3*NRATMS, TOTREDNCO, 1)
#endif
      CALL DGETREC(20,'JOBARC','GMATRIX',LENGTH_BGMAT,
     &            GMINBT)

#ifdef _DEBUG_LVL0
         print *, 'The G matrix'
         CALL OUTPUT(GMINBT, 1, 3*NRATMS, 1, TOTREDNCO, 
     &               3*NRATMS, TOTREDNCO, 1)
#endif
C   Here we transpose the GMATRIX from JOBARC, because it was stored as
C   (G(-)*B)t and we need (G(-)*B)

      CALL TRANSP(GMINBT,GMINB,TOTREDNCO,3*NRATMS)

C   And here what we do is G(-)*B*(HC-K)*B(t)*G(-)

      CALL XGEMM('N','N',TOTREDNCO,3*NRATMS,3*NRATMS,1.0D0,
     &           GMINB,TOTREDNCO,DIFTEMP,3*NRATMS,0.0D0,
     &           TEMP1,TOTREDNCO)

      CALL XGEMM('N','N',TOTREDNCO,TOTREDNCO,3*NRATMS,1.0D0,
     &           TEMP1,TOTREDNCO,dRICHeap(z_BTGInv),3*NRATMS,0.0D0,
     &           HI,TOTREDNCO)

C Now we have the the Hessian in Redundant Internal Coordinates
#ifdef _DEBUG_LVL0
         print *, 'Hessian in redundant internals'
         CALL OUTPUT(HI, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
     &               TOTREDNCO, 1)
#endif
C Purify the Hessian as recommended by Schlegel et al. (JCC,17,49,1996).
C G and G{^-1} have been created in the built_bgmatrix.F. The P=GG{^-1}
C is the projection into the redundant spaces that correspond to actual
C Cartesian displacements.

      LENGMAT=TOTREDNCO*TOTREDNCO
      CALL DGETREC(20,'JOBARC','G-MATRX ',LENGMAT,DIFTEMP)
      CALL DGETREC(20,'JOBARC','GI-MATRX',LENGMAT,
     &            dRICHeap(z_BTGInv))

      CALL ZERO(TEMP1,LENGMAT)
      CALL XGEMM('N','N',TOTREDNCO,TOTREDNCO,TOTREDNCO,1.0D0,
     &           DIFTEMP,TOTREDNCO,dRICHeap(z_BTGInv),TOTREDNCO,0.0D0,
     &           TEMP1,TOTREDNCO)

      CALL DPUTREC(20,'JOBARC',"PROJECTR",LENGMAT,TEMP1)

      CALL XGEMM('N','N',TOTREDNCO,TOTREDNCO,TOTREDNCO,1.0D0,
     &           TEMP1,TOTREDNCO,HI,TOTREDNCO,0.0D0,GMINB,
     &           TOTREDNCO)
      CALL XGEMM('N','N',TOTREDNCO,TOTREDNCO,TOTREDNCO,1.0D0,
     &           GMINB,TOTREDNCO,TEMP1,TOTREDNCO,0.0D0, HI,
     &           TOTREDNCO)

C Also project the gradients. Then copy the projected gradient
C to the FI array so that the rest of the code can proceed.
      CALL XGEMM('N', 'N',TOTREDNCO,1,TOTREDNCO,1.0D0,
     &           TEMP1,TOTREDNCO,FI,TOTREDNCO,0.0D0,
     &           PFI,TOTREDNCO)
      CALL XDCOPY(TOTREDNCO,PFI,1,FI,1)

C Use GMINB to store the identity matrix.
      CALL ZERO(GMINB,TOTREDNCO*TOTREDNCO)
      DO I=1,TOTREDNCO
         GMINB((I-1)*TOTREDNCO+I) = 1.d0
      END DO

C Let's do the PHP + 1000.0(1-P) (note that we have already
C built the PHP and stored in HI).
      CALL XDAXPY(TOTREDNCO*TOTREDNCO,-1.d0,GMINB,1,TEMP1,1)

      CALL XDSCAL(TOTREDNCO*TOTREDNCO,1000.0D0,TEMP1,1)
      CALL XDAXPY(TOTREDNCO*TOTREDNCO,-1.d0,TEMP1,1,HI,1)
#ifdef _DEBUG_LVL0
          print*,'The projected hessian in redundent internals '
          CALL OUTPUT(HI,1,TOTREDNCO,1,TOTREDNCO,TOTREDNCO,
     &                TOTREDNCO,1)
#endif
      RETURN
      END

