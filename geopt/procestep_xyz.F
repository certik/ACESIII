      SUBROUTINE PROCESTEP_XYZ(SCRATCH, AMATRX, TOTREDNCO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL CONVERGED, CHANGED, NO_ITER
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG
      DOUBLE PRECISION ORIENT(3,3)
      PARAMETER (MAXITER = 100)
C
#include "mxatms.par"
#include "machsp.com"
#include "flags.h"
#include "flags.com"
#include "flags2.com"
#include "coord.com"
C
      PARAMETER (EPSILON = 1.0D-5)
      CHARACTER*80 HOLDQ
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION SCRATCH(9*NATOMS*NATOMS), IQFIX(3*MXATMS, 3*MXATMS),
     &          AMATRX(3*NATOMS*TOTREDNCO), STATS(6), 
     &          REDUNCO(MAXREDUNCO), SCRATCH_LCL(4*MAXREDUNCO+6*
     &          MXATMS+15*MXATMS),
     &          TMP1(MAXREDUNCO*MAXREDUNCO),TMP2(MAXREDUNCO)
C
C TMP1(17*TOTREDNCO),TMp2(17)
C
      DATA BTOA /0.529177249D0/

      DO I = 1, NUNIQUE
         DO J = 1, NEQ(IUNIQUE(I))
            IQFIX(IUNIQUE(I), J) = IEQUIV(I,J)
         ENDDO
      ENDDO
C
      DO J = 1, NOPT
C
         SCRATCH_LCL(NOPTI(J)) = SCRATCH(J + NOPT)
         SCRATCH_LCL(NOPT + NOPTI(J)) = SCRATCH(J + NOPT)
C
         DO K = 1, NEQ(NOPTI(J))
            SCRATCH_LCL(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
           SCRATCH_LCL(NOPT + IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)  
         ENDDO
      ENDDO
C
#ifdef _DEBUG_LVLM2
      Write(6,*)
      WRITE(6,*) "The optimization parameters:nx, nopt, totrednco nxm6",
     &            nx, totrednco, nopt, nxm6
      Write(6,*)
      WRITE(6,*) "The delta q0"
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(I), I = 1, TOTREDNCO)
   10 Format (5(1X,F10.6))
      Write(6,*)
#endif
C
      IOFF4q0  = 3*NATOMS + 3*TOTREDNCO
      IOFF4X0  = IOFF4q0  + TOTREDNCO
C
      ISCRATCH = 4*TOTREDNCO + 6*NATOMS 
      ITMPQ    = ISCRATCH    + 9*NATOMS
      INEWQ    = ITMPQ       + 3*NATOMS
      ILAST    = INEWQ       + 3*NATOMS
C
      CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL, 1, R, 1)
      CALL DCOPY(TOTREDNCO, R, 1, SCRATCH_LCL(2*TOTREDNCO + 1), 1)
      CALL DCOPY(TOTREDNCO, SCRATCH_LCL, 1, SCRATCH_LCL(IOFF4q0 + 1),
     &           1)
      CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(IOFF4X0 + 1), 1)

C#ifdef _TEST_RED
C This requires defining TMP1 and TMP2 arrays.
      CALL IGETREC(20,'JOBARC','NUMREDCO', 1, NULLEVAL)
      CALL DGETREC(20,'JOBARC','REDEVECS',NULLEVAL*TOTREDNCO,
     &             TMP1)
      CALL XGEMM('T', 'N', NULLEVAL, 1, TOTREDNCO, 1.0D0,
     &               TMP1, TOTREDNCO, SCRATCH_LCL, TOTREDNCO,
     &               0.0D0, TMP2, NULLEVAL)
#ifdef _DEBUG_LVLM1
      Write(6,*) "Check for the redundencies in deltaq0"
      Write(6,"(6F13.7)") (tmp2(i), i=1, NULLEVAL)
#endif
C#endif
C
#ifdef _DEBUG_LVLM1
      WRITE(6,*) "The starting Cartesian coords. x0"
      Write(6,*)
      CALL OUTPUT(Q, 1, 3, 1, Natoms, 3, Natoms, 1)
      Write(6,*) 
      Write(6,*) "The starting RICs  (q0)"
      Write(6,*) 
      WRITE(6,10) (SCRATCH_LCL(2*TOTREDNCO + I), I=1, TOTREDNCO)
   10 Format (5(1X,F10.6))
      Write(6,*)
#endif
C     
      CONVERGED = .FALSE. 
        CHANGED = .FALSE.
        NO_ITER = .FALSE.
          NITER = 0 
C
      DO WHILE (.NOT. CONVERGED .AND. NITER .LE. MAXITER .AND.
     &          .NOT. CHANGED .AND. .NOT. NO_ITER) 
C
         CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(3*TOTREDNCO + 1), 1)
#ifdef _DEBUG_LVLM1
         Write(6,*)
CSSS         Write(6,*) "The A matrix"
CSSS         CALL OUTPUT(AMATRX, 1, 3*NATOMS, 1, TOTREDNCO, 3*NATOMS,
CSSS     &               TOTREDNCO, 1) 
         WRITE(6,*) "The starting Cartesian during iteration; x_K"
         CALL OUTPUT(Q, 1, 3, 1, Natoms,3, Natoms, 1)
         Write(6,*)
         Write(6,*) "The starting internal increments:delta q_k"
         Write(6,*) 
         WRITE(6,10) (SCRATCH_LCL(I), I=1, TOTREDNCO)
         Write(6,*)
#endif
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &               AMATRX, 3*NATOMS, SCRATCH_LCL, TOTREDNCO,
     &               1.0D0, Q, 3*NATOMS)
C
CSSS         CALL DGETREC(20,'JOBARC','ORIENT2 ',9*IINTFP,ORIENT)
CSSS         CALL XGEMM('T', 'N', 3, NATOMS, 3, 1.0D0, ORIENT, 3, Q,
CSSS     &               3, 0.0D0, SCRATCH_LCL(INEWQ), 3)
CSSS         CALL  DCOPY(3*NATOMS, SCRATCH_LCL(INEWQ), 1, Q, 1)
C
#ifdef _DEBUG_LVL0
         Print*, "The Cartesian after the update: x_k+1"
         CALL OUTPUT(Q, 1, 3, 1, Natoms,3, Natoms, 1)
#endif
         CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNOFBND,
     &                    NW_TOTREDNCO, TOTNOFANG, .FALSE.)
C
         IF (IFLAGS(h_IFLAGS_symmetry).EQ. 2) THEN
            CALL DGETREC(20,'JOBARC','PLSMINSP',NW_TOTREDNCO,NCON)  
            DO ICOORD = 1, NW_TOTREDNCO
               IF (NCON(ICOORD).EQ.1) REDUNCO(ICOORD) = -REDUNCO(ICOORD)
            END DO
         END IF
C
C In rare cases, it is possible that the # of RIC coordinates
C can change. If that happens leave the iterative procedure and
C and use the linear update.
C
         IF (NW_TOTREDNCO .GT. TOTREDNCO) CHANGED = .TRUE.
         CALL CHECK4_180(SCRATCH_LCL(2*TOTREDNCO + 1), SCRATCH_LCL,
     &                   TOTNOFBND, TOTNOFANG, TOTREDNCO, NO_ITER)
C     
#ifdef _DEBUG_LVLM1
      Write(6,*)
      WRITE(6,*) "The new q: q_k+1"
      Write(6,*)
      WRITE(6,10) (REDUNCO(I), I = 1, TOTREDNCO)
      Write(6,*) 
      WRITE(6,*) "The strating q: q0 (should not change)"
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(2*TOTREDNCO+I), I = 1, TOTREDNCO)
      Write(6,*)
#endif
C
         CALL DCOPY(TOTREDNCO, REDUNCO, 1, R, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL(2*TOTREDNCO + 1),
     &              1, R, 1)
#ifdef _DEBUG_LVLM1 
      WRITE(6,*) "The (q_k - q0)"
      Write(6,*)
      WRITE(6,10) (R(I), I = 1, TOTREDNCO)
      Write(6,*)
      Write(6,*) "The delta q0 (should not change)" 
      Write(6,*)
      write(6,10) (SCRATCH_LCL(IOFF4q0 + I), I=1, TOTREDNCO)
      Write(6,*)
#endif
         CALL DCOPY(TOTREDNCO, SCRATCH_LCL(IOFF4q0 + 1), 1, 
     &              SCRATCH_LCL, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, R, 1, SCRATCH_LCL, 1)
C
#ifdef _DEBUG_LVLM1
      WRITE(6,*) "The delta delta q_k" 
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(I), I=1, TOTREDNCO)
      Write(6,*)
#endif
         CALL REMOVE_360(SCRATCH_LCL, TOTNOFBND, TOTNOFANG, 
     &                   TOTREDNCO)
C
         DO IREDUNCO = 1, TOTREDNCO
            IF (DABS(SCRATCH_LCL(IREDUNCO)) .LT. 1.0D-10) SCRATCH_LCL
     &          (IREDUNCO) = 0.0D0
         ENDDo
   
         IF (NITER .GE. 0) CALL VSTAT(SCRATCH_LCL, STATS, TOTREDNCO)
         IF (NITER .EQ. 0) THEN
             RMSOF_DDQ1 = STATS(5)
         ELSE
C
#ifdef _DEBUG_LVLM1
      Write(6,"(a,2F13.9)"), "RMS RIC changes DDQ1,DDDQk:", 
     &     RMSOF_DDQ1, STATS(5)
#endif
CSSS             IF (STATS(5) .GE. RMSOF_DDQ1) NO_ITER = .TRUE.
             IF (STATS(5) .GE. RMSOF_DDQ1) NO_ITER = .FALSE.
         ENDIF
C
         CALL DAXPY(3*NATOMS, -1.0D0, Q, 1, SCRATCH_LCL
     &             (3*TOTREDNCO + 1), 1)
         CALL VSTAT(SCRATCH_LCL(3*TOTREDNCO+1), STATS, 3*NATOMS)
C
#ifdef _DEBUG_LVLM1
         Write(6,*)
         Write(6,*) "Change in Cartesians x_k+1 - x_k"
         CALL OUTPUT(SCRATCH_LCL(3*TOTREDNCO+1), 1, 3, 1, NATOMS,
     &               3, NATOMS, 1)
         Write(6,*) 
         Print*, "The root mean sq change of x:", STATS(5)
#endif
         IF (STATS(5)  .LE. EPSILON) CONVERGED = .TRUE.
         NITER = NITER + 1 

      ENDDO 
      IF (CONVERGED) Write(6, "(a,a,I2,a)") " RIC iterative update ",
     &                         "converged in ", NITER, " iterations."
C
C If the iterative procedure fail to converge in maximum iterations
C alowed, then do a simple linear update.
C
C#endif
C
      IF (.NOT. CONVERGED) THEN
          Write(6, "(a,a)") " RIC iterative update did not converged.",
     &                      " linear update is used. "
         CALL DCOPY(3*NATOMS, SCRATCH_LCL(IOFF4X0 +1), 1, Q, 1)
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &              AMATRX, 3*NATOMS, SCRATCH_LCL(IOFF4q0 + 1), 
     &              TOTREDNCO, 1.0D0, Q, 3*NATOMS)
      ENDIF
C
C Lets screen the Cartesian coordiantes and filter out changes 
C beyond 7th decimal. 
C
      DO IATOM = 1, 3*NATOMS
         WRITE(HOLDQ, "(F12.7)") Q(IATOM)
         READ(HOLDQ,  "(F12.7)") Q(IATOM)
      ENDDO
C
#ifdef _DEBUG_LVLM1
      WRITE(6,*) "The new Cartesian Coord., x+delta x0"
      WRITE(6,"(3F10.5)") (Q(I), I = 1, 3*NATOMS)
#endif

C
C Store the Cartesian coordinates with a different record label
C from COORD since finite difference calculations use COORD
C for the grid points, not the current point.
C
      CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNOFBND,
     &                 NW_TOTREDNCO, TOTNOFANG, .TRUE.)

      CALL DPUTREC(1,'JOBARC','COORD_OP',3*NATOMS,Q)
C
      RETURN
      END
