      SUBROUTINE FLOWEVFMIN(SCRATCH, GRDHES, HESMOD, DIAGHES, LMBDAN,
     &                      MORSE, NX, NOPT)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DOUBLE PRECISION LMBDAN
      LOGICAL MORSE
C
#include "flags.h"
#include "jodaflags.com"
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), HESMOD(NOPT, NOPT),
     &          DIAGHES(NOPT, NOPT)
C
#ifdef _DEBUG_LVLM1
        Write(6,*) "The scratch array before RFO step", NOPT
        Write(6, "(3F10.5)") (SCRATCH(J+NOPT), J=1, NOPT)
#endif
      DO 10 I = 1, NOPT
         DO 20 J = 1, NOPT

            IF (LAMBDN .EQ. 0.0D0 .AND.
     &          DABS(HESMOD(I,I)) .LT. 1.0D0-8)  THEN
C
                 SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
C
            ELSE
C
               IF (DABS((HESMOD(I,I) - LMBDAN)) .LT. 1.0D0-8) THEN

                  SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
C
               ELSE
                  SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-
     &                             GRDHES(I)*DIAGHES(J,I)/
     &                             (HESMOD(I,I) - LMBDAN)
C
#ifdef _DEBUG_LVLM1
            Write(6,"(6F10.5)"), GRDHES(I), HESMOD(i,i),
     &                    DIAGHES(J,I), LMBDAN, GRDHES(I)/
     &                    (HESMOD(I,I) - LMBDAN), SCRATCH(J+NOPT)
#endif
              ENDIF
            ENDIF

 20      CONTINUE
 10   CONTINUE 
#ifdef _DEBUG_LVLM1
        Write(6,*) "The unscaled RFO step"
        Write(6, "(3F10.5)") (SCRATCH(J+NOPT), J=1, NOPT)
#endif
C
C Notice that there is no Morse adjustments can be made for
C pure Cartesian OLPtimization (no connections), Also, no
C Morse scaling can be done for redundent internals.
C Ajith Perera,08/2008, 2011.

      IF (MORSE) THEN
         IF (iFlags2(h_IFLAGS2_geom_opt) .ge. 3) THEN
CSSS           CALL DOMORSEXYZ_RIC(SCRATCH, NOPT, NX, LUOUT)
         ELSE IF (iFlags2(h_IFLAGS2_geom_opt) .eq. 1) THEN
           CALL DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
         ENDIF
      ENDIF
C
      RETURN
      END
