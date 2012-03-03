#include "flags.h"
      SUBROUTINE A2MAKE_WFN(DENSA, DENSB, NATORBA, NATORBB, EVECA,
     &                      EVECB, EIGVALA, EIGVALB, SCR1, SCR2,
     &                      OCCA, OCCB, NOCCA, NOCCB, MAXOCCA,
     &                      MAXOCCB, NAOBFNS, NBAS, IUHF, TENERGY)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "machsp.com"
#include "symm2.com"
#include "jodaflags.com"
C
      DOUBLE PRECISION NATORBA, NATORBB, OCCA(NBAS), OCCB(NBAS)
      LOGICAL SCF, COR
      DIMENSION DENSA(NBAS*NBAS), DENSB(NBAS*NBAS), NOCC(16),
     &          NATORBA(NAOBFNS, NAOBFNS), NATORBB(NAOBFNS, NAOBFNS),
     &          EVECA(NBAS, NBAS), EVECB(NBAS, NBAS),
     &          EIGVALA(NBAS), EIGVALB(NBAS), 
     &          SCR1(NBAS*NBAS), SCR2(NBAS*NBAS)
C
      DATA ZILCH, ONE, HALF/ 0.0D0, 1.0D0, 0.50D0/
C
      SCF = (iflags(h_IFLAGS_calc) .eq.0)
      COR = (iflags(h_IFLAGS_calc) .gt.0)
C
      
      IF (SCF) THEN

         CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP, EVECA)
         CALL GETREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP, EVALA)
         DO IOCC = 1, NOCCA
            OCCA(IOCC) = 2.0D0
         ENDDO
         MAXOCCA = NOCCA
C
         CALL GETREC(20, "JOBARC", "KINETINT", NBAS*NBAS*IINTFP, SCR1)
         CALL GETREC(20, "JOBARC", "SCFDENSA", NBAS*NBAS*IINTFP, SCR2)
         TENERGYA = DDOT(NBAS*NBAS, SCR1, 1, SCR2, 1)
C
         IF (IUHF .EQ. 1) THEN
            CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,
     &                  EVECB)
            CALL GETREC(20,'JOBARC','SCFEVLB0',NBAS*IINTFP, EVALB)
            Do IOCC = 1, NOCCA
               OCCA(IOCC) = 1.0D0
            Enddo
            Do IOCC = 1, NOCCB
                OCCB(IOCC) = 1.0D0
            ENDDO
            MAXOCCB = NOCCB
C
         CALL GETREC(20, "JOBARC", "SCFDENSB", NBAS*NBAS*IINTFP, SCR2)
         TENERGYB = DDOT(NBAS*NBAS, SCR1, 1, SCR2, 1)
C
         ENDIF
         TENERGY = TENERGYA + TENERGYB
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,F10.5)") "The total kinetic energy = ", TENERGY
#endif 
C
      ELSE IF (COR) THEN
C
         CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,DENSA)
         CALL DCOPY(NBAS*NBAS, DENSA, 1, SCR2, 1)
         CALL MO2AO3(SCR2, DENSA, EVECA, SCR1, NBAS, NBAS, 1, 2)
         CALL GETREC(20, "JOBARC", "KINETINT", NBAS*NBAS*IINTFP, SCR1)
         TENERGYA = DDOT(NBAS*NBAS, DENSA, 1, SCR1, 1)
C
         CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,DENSA)
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "The Alpha Density read"
      call output(DENSA, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
C
         CALL EIG(DENSA, SCR1, 1, NBAS, -1)
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "The Alpha Natural orbitals MOxMO"
      call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
         CALL DCOPY(NBAS, DENSA, NBAS+1, OCCA, 1)
         CALL DZERO(EVALA, NBAS)
         MAXOCCA = NBAS
         CALL GETREC(20,'JOBARC','SCFEVECA',NBAS*NBAS*IINTFP,DENSA)
         CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,DENSA,NBAS,SCR1, 
     &               NBAS,ZILCH,EVECA,NBAS)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The Alpha Natural orbitals AOxMO"
      call output(EVECA, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
C
        IF (IUHF.EQ.1) THEN
           CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,DENSB)
           CALL DCOPY(NBAS*NBAS, DENSB, 1, SCR2, 1)
           CALL MO2AO3(SCR2, DENSB, EVECA, SCR1, NBAS, NBAS, 2, 2)
           CALL GETREC(20, "JOBARC", "KINETINT", NBAS*NBAS*IINTFP, SCR1)
           TENERGYB = DDOT(NBAS*NBAS, DENSB, 1, SCR1, 1)
C
           CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,DENSB)
C
#ifdef _DEBUG_LVL0
       Write(6,*) "The Beta Density read"
       call output(DENSB, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
           CALL EIG(DENSB, SCR1, 1, NBAS, -1)
           CALL DCOPY(NBAS, DENSB, NBAS+1, OCCB, 1)
           CALL DZERO(EVALB, NBAS)
           MAXOCCB = NBAS
           CALL GETREC(20,'JOBARC','SCFEVECB',NBAS*NBAS*IINTFP,DENSB)
           CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,DENSB,NBAS,SCR1,
     &                 NBAS,ZILCH,EVECB,NBAS)
#ifdef _DEBUG_LVL0
      Write(6,*) "The Beta Natural orbitals AOxMO"
      call output(EVECB, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
        ENDIF
        TENERGY = TENERGYA + TENERGYB

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,F10.5)") "The total kinetic energy = ", TENERGY
#endif 
C
      ENDIF
C
      CALL GETREC(20,'JOBARC','CMP2ZMAT', NBAS*NAOBFNS*IINTFP, SCR1)
      CALl XGEMM('N','N', NAOBFNS, NBAS, NBAS, ONE, SCR1, NAOBFNS, 
     &            EVECA, NBAS, ZILCH, NATORBA, NAOBFNS)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The Alpha Natural orbitals NAOBFNSxNMO"
      call output(NATORBA, 1, NAOBFNS, 1, NBAS, NAOBFNS,  NBAS, 1)
#endif
C
      IF (IUHF .EQ. 1) THEN
         CALL XGEMM('N','N', NAOBFNS, NBAS, NBAS, ONE, SCR1, NAOBFNS,
     &               EVECB, NBAS, ZILCH, NATORBB, NAOBFNS)
      ENDIF
  
      RETURN
      END
