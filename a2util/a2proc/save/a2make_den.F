
      SUBROUTINE A2GET_DEN(TMPDEN, TOTDEN, SPNDEN)
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPHERICAL, SPIN_D
C
      CHARACTER*8 LABELPSCF, LABELDENS, LABELD,
      DIMENSION TMPDEN(NBAS*NBAS), TOTDEN(NBAS*NBAS), 
     &          SPNDEN(NBAS*NBAS)
C
      DATA ONE, ONEM /1.0D0, -1.0D0/
C
#include "machsp.com"
#include "syminf.com"
#include "jodaflags.com"
C
      LABELPSCF ='HFDENSTY'
      LABELDENS ='TDENSITY'
C
C Get the SCF and correlated (relaxed) density matirices in AO basis. For 
C open shell calculations, we also read the density diferences. Note that
C the "relaxed" density read here only contain the correlated contributions.

      IF (SCF) THEN

         DO ISPIN = 1, (IUHF + 1)
C            
           IF (ISPIN .EQ. 1) THEN
               LABELD='SCFDENSA'
            ELSE
               LABELD='SCFDENSB'
            ENDIF
C           
            CALL GETREC(20, 'JOBARC', LABELD, NBAS*NBAS*IINTFP, TMPDEN)
            IF (ISPIN .EQ. 1) THEN
               CALL DCOPY(NBAS*NBAS, TMPDEN, 1, TOTDEN, 1)
            ELSE
               CALL SAXPY(NBAS*NBAS, ONE, TMPDEN, 1, TOTDEN, 1)
            ENDIF
C
            IF (IUHF .NE. 0) THEN
               IF (ISPIN .EQ. 1) THEN
                  CALL DCOPY(NBAS*NBAS, TMPDEN, 1, SPNDEN, 1)
               ELSE
                  CALL SAXPY(NBAS*NBAS, ONEM, EVEC, 1, SPNDEN, 1)
               ENDIF
            ENDIF
C
         ENDDO
C
C        
      ELSE IF (COR) THEN

         DO ISPIN = 1, (IUHF + 1)
C   
           IF (ISPIN .EQ. 1) THEN
               LABELD='RELDENSA'
            ELSE
               LABELD='RELDENSB'
            ENDIF
C
            CALL GETREC(20, 'JOBARC', LABELD, NBAS*NBAS*IINTFP, TMPDEN)
            IF (ISPIN .EQ. 1) THEN
               CALL DCOPY(NBAS*NBAS, TMPDEN, 1, TOTDEN, 1)
            ELSE
               CALL SAXPY(NBAS*NBAS, ONE, TMPDEN, 1, TOTDEN, 1)
            ENDIF
C
            IF (IUHF .NE. 0) THEN
               IF (ISPIN .EQ. 1) THEN
                  CALL DCOPY(NBAS*NBAS, TMPDEN, 1, SPNDEN, 1)
               ELSE
                  CALL SAXPY(NBAS*NBAS, ONEM, EVEC, 1, SPNDEN, 1)
               ENDIF
            ENDIF
         ENDDO

      ENDIF
C
      RETURN
      END


