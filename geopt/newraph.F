      SUBROUTINE NEWRAPH(SCRATCH, GRDHES, HESMOD, DIAGHES, MORSE, 
     &                   NOPT, NX, LUOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "coord.com"
#include "flags.h"
#include "jodaflags.com"
C
      LOGICAL MORSE, XYZIN, NWFINDIF
C
      COMMON /INPTYP/ XYZIN,NWFINDIF
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), DIAGHES(NOPT, NOPT),
     &          HESMOD(NOPT, NOPT)
C     
C Do the Normal Newton-Raphson update
C
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6, "(1x,a,I5)") "The number of degs. of freed. 
     &          at NR:", NOPT
      Write(6,*)
#endif
      DO 10 I = 1, NOPT
         DO 20 J = 1, NOPT
            IF (HESMOD(I,I) .EQ. 0.0D0) Then
                    SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
            ELSE
C
                SCRATCH(J+NOPT)=SCRATCH(J+NOPT)-GRDHES(I)*
     &                          DIAGHES(J,I)/HESMOD(I,I) 
            ENDIF
 20      CONTINUE
 10   CONTINUE
C
C Notice that there is no Morse adjustments can be made for
C pure Cartesian OLPtimization (no connections), Also, no
C Morse scaling can be done for redundent internals.
C Ajith Perera,08/2008, 2011.
C
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
