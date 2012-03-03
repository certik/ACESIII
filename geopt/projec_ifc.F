      Subroutine projec_IFC(HESS, PHESS, PMAT, GMATRX_N, GMATRX_M,
     &                      GRD, NXM6)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include <machsp.com>
C
      Dimension HESS(NXM6, NXM6), PMAT(NXM6, NXM6), PHESS(NXM6, NXM6),
     &          GMATRX_M(NXM6, NXM6), GMATRX_N(NXM6, NXM6),
     &          GRD(NXM6)
 
      LENGMAT=NXM6*NXM6
      CALL DGETREC(20,'JOBARC','GI-MATRX',LENGMAT,GMATRX_M)
      CALL DGETREC(20,'JOBARC','G-MATRX ',LENGMAT,GMATRX_N)

      Call XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, GMATRX_M, 
     &            NXM6, GMATRX_N, NXM6, 0.0D0, PMAT, NXM6)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The Hessian projector:P"
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)

#endif
C
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)

#ifdef _DEBUG_HESSIAN_PROJC
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The (PHP)(PHP)t matrix"
      CALL OUTPUT(HESS, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
#endif
C
      CALL EIG(HESS, PHESS, NXM6, NXM6, 1)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,"(a)") "The eigen vectors of the PHP"
      CALL OUTPUT(phess, 1, nxm6, 1, nxM6, nxm6, nxm6, 1)
      Write(6,"(a)") "The eigenvalues of the PHP "
      Write(6, "(4F10.5)") (Hess(I,I), I=1, Nxm6)
#endif
      DO I = 1, NXM6
         IF (HESS(I, I) .LE. 1.0D-10) THEN
             HESS(I, I) = 0.0D0
         ELSE
             HESS(I, I) = 1.0D0/HESS(I, I)
         ENDIF
      ENDDO
C
C Built the generalized inverse of H-matrix,
C
      CALL XGEMM('N', 'N', NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, HESS, NXM6, 0.0D0, GMATRX_N, NXM6)
      CALL XGEMM('N', 'T', NXM6, NXM6, NXM6, 1.0D0, GMATRX_N,
     &            NXM6, PHESS, NXM6, 0.0D0, HESS, NXM6)
C
C Built P(PHP)-P 
C
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)
      CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)
C
C---These steps require defining STEP and TMP arrays.
C
      CALL XGEMM("N", "N", NXM6, 1, NXM6, -1.0D0, HESS,
     &            NXM6, GRD, NXM6, 0.0D0, STEP, NXM6)
      NULLEVAL=8
      CALL DGETREC(20,'JOBARC','REDEVECS',NULLEVAL*TOTREDNCO,
     &            GMATRX_N)
      CALL XGEMM('T', 'N', 8, 1, TOTREDNCO, 1.0D0,
     &               GMATRX_N, TOTREDNCO, R, TOTREDNCO,
     &               0.0D0, TMP, 8)
     
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The step and overlap with the red. space"
      Write(6, "(8F10.5)") (Step(I), I=1, Nxm6)
      Write(6, "(8F10.5)") (TMP(I), I=1, 8)
#endif
      STOP
#endif
      RETURN
      END
     
