
C Compute a "new" gradient using the SCF gradient and the
C total (SCF + correlation) gradient. This can be considered
C as one (among a family) of gradient extrapolation schemes.
C See (???) for further details.

      SUBROUTINE EXTRAP_GRAD(SCF_GRAD, TOT_GRAD, NATOMS, STRING)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "machsp.com"
#include "flags.h"
#include "flags.com"
C
      DIMENSION SCF_GRAD(3*NATOMS), TOT_GRAD(3*NATOMS)
      LOGICAL bExist, bPrint
      CHARACTER*(*) STRING
C
      bPrint = (iflags(h_IFLAGS_print).GE.10)
C
C Takashi, extrapolation formula goes here, and the extrapolated
C gradient should be stored back into TOT_GRAD.
C
      INQUIRE(FILE='fact',EXIST=bExist)
      IF (bExist) THEN
c         OPEN(99,FILE='fact',FORM='UNFORMATTED', STATUS='OLD')
         OPEN(UNIT=99,FILE='fact',FORM='FORMATTED',STATUS='OLD')
         READ(99,*) FACT
         CLOSE(99)
c         write(*,*) FACT
      ELSE
         PRINT *, "@EXTRAP_GRAD: The formatted file 'fact' is missing."
         CALL ERREX
      END IF
C
      IF (STRING.EQ.'gradient') THEN
C
         CALL GETREC(1, "JOBARC", "SCF_GRAD", 3*NATOMS*IINTFP, SCF_GRAD)
         CALL GETREC(1, "JOBARC", "GRADIENT", 3*NATOMS*IINTFP, TOT_GRAD)
         IF (bPrint) THEN
            write(*,'(" SCF_GRAD:\n",3(F15.10))')
     &               (SCF_GRAD(J),J=1,3*NATOMS)
            write(*,'(" TOT_GRAD:\n",3(F15.10))')
     &               (TOT_GRAD(J),J=1,3*NATOMS)
         END IF
C
         DO I = 1,3*NATOMS
            TOT_GRAD(I) =   SCF_GRAD(I)
     &                    + (TOT_GRAD(I)-SCF_GRAD(I))/FACT
         END DO
         IF (bPrint) THEN
            PRINT *
            PRINT *, " The Extrapolated Total Gradients"
            CALL OUTPUT(TOT_GRAD, 1, NATOMS, 1, 3, NATOMS, 3, 0)
         END IF
C
         CALL PUTREC(1, "JOBARC", "GRADIENT", 3*NATOMS*IINTFP, TOT_GRAD)

      ELSE IF (STRING.EQ."energy") THEN
C
         CALL GETREC(1, "JOBARC", "SCFENEG ", IINTFP, SCF_ENEG)
         CALL GETREC(1, "JOBARC", "TOTENERG", IINTFP, ETOT)
         IF (bPrint) THEN
            write(*,'(" SCF_ENEG=",F16.10," a.u.")') SCF_ENEG
            write(*,'(" ETOT    =",F16.10," a.u.")') ETOT
         END IF
C
         ESAC = SCF_ENEG + (ETOT-SCF_ENEG)/FACT
         WRITE(6,2) ESAC
    2    FORMAT('  ','Total MBPT-SAC(2) energy =',F18.12,' a.u.')
C
         CALL PUTREC(1, "JOBARC", "TOTENERG", IINTFP, ESAC)
C
      ELSE

         print *, '@EXTRAP_GRAD: unknown argument'
         print *, '              TYPE = ',string(1:20)
         call errex

      END IF

      RETURN
      END

