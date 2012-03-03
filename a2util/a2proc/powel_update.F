
      SUBROUTINE POWEL_UPDATE(VC,VP,H,SCRATCH,STEP,TBT,PHI,NXM6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION VC(NXM6),H(NXM6*NXM6),SCRATCH(NXM6*NXM6),
     &          TBT(3*NXM6*NXM6),STEP(NXM6),VP(NXM6)
C
C X = 1/DX{dot}DX where DX is the step size.
C
C      X=1.D0/xdot(NXM6,STEP,1,STEP,1)
      X=xdot(NXM6,STEP,1,STEP,1)
      If ((1.0D0/X) .LT. 1.0D-08) Return

      N2=1
      N3=NXM6+1
      N4=2*NXM6+1
C
C SCRATCH(N2) = HDX where H is the Hessian, SCRATCH(N3) = DG (gradient
C difference), SCRATCH(N4) = (DG - HDX)
C
      CALL MODMATMUL(SCRATCH(N2),H,STEP,NXM6,NXM6,1,NXM6,NXM6,1)
      CALL VADD(SCRATCH(N3),VC,VP,NXM6,-1.D0)
      CALL VADD(SCRATCH(N4),SCRATCH(N3),SCRATCH(N2),NXM6,-1.D0)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "@-Powel_update DG and HDX and DG-HDX"
      Write(6, "(3F12.8)") (scratch((n3-1)+i),i=1,NXM6)
      Write(6,*)
      Write(6, "(3F12.8)") (scratch((n2-1)+i),i=1,NXM6)
      Write(6,*)
      Write(6, "(3F12.8)") (scratch((n4-1)+i),i=1,NXM6)
#endif
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
     
      CALL xscal(NXM6*NXM6,dtmp/X,TBT(1),1)
#ifdef _DEBUG_LVL0
      Write(6,"(a,3F17.13)")"The scaling factor",dtmp,dtmp/X, x
#endif
C
C TBT(N3) = (DG - HDX)DX + [(DG - HDX)DX]^{t} - [DXdot(DG - HDX)]DXDX/DX{dot}DX
C
      CALL VADD(TBT(N3),TBT(N3),TBT(N2),NXM6*NXM6,1.D0)
      CALL VADD(TBT(N3),TBT(N3),TBT(1),NXM6*NXM6,-1.D0)
C
C H(updated) = H(old) + [(DG - HDX)DX + [(DG - HDX)DX]^{t}]/DX{dot}DX +
C              DXdot(DG - HDX)]DXDX/[DX{dot}DX]^{2}
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "Intermediate check"
      call output(TBT(N3), 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif

      CALL xscal(NXM6*NXM6,1.0D0/x,TBT(N3),1)
      CALL DAXPY(NXM6*NXM6, PHI, TBT(N3), 1, H, 1)
C
C DON'T WRITE ENTIRE HESSIAN UNLESS SPECIFICALLY REQUESTED.  FOR
C NOW, JUST USE AN IN-CODE PARAMETER.
C
#ifdef _DEBUG_LVL0
         Print*, "-----The Powel updated Hessian-----"
         CALL output(H,1,NXM6,1,NXM6,NXM6,NXm6,1)
#endif
      RETURN
      END

