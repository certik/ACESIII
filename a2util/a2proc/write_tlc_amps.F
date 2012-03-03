      SUBROUTINE WRITE_TLC_AMPS(T, TLEN, NAME)
C
C This routine dumps the T vector into a file to be picked up and used
C as a guess in a later CC calculation.
C
C T contains the T amplitudes in the following order
C    T1AI   (List 90,1)
C    T1ai   (List 90,2) (if IUHF .NE. 0)
C    T2abij (List 45)   (if IUHF .NE. 0)
C    T2ABIJ (List 44)
C    T2AbIj (List 46)
C
      IMPLICIT NONE
C
      INTEGER TLEN
      DOUBLE PRECISION T(TLEN)
      CHARACTER *8 NAME
C
      INTEGER NAMLEN
      LOGICAL YESNO
      CHARACTER *80 FULNAM
C
      CALL GFNAME(NAME, FULNAM, NAMLEN)
C
      INQUIRE(FILE=FULNAM(1:NAMLEN), EXIST=YESNO)
      IF (YESNO) THEN
        OPEN(UNIT=94, FILE=FULNAM(1:NAMLEN), STATUS='OLD',
     &     FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
        CLOSE(UNIT=94, STATUS='DELETE')
      ENDIF
C
      OPEN(UNIT=94, FILE=FULNAM(1:NAMLEN), STATUS='NEW',
     &   FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
      WRITE(94) T
      CLOSE(UNIT=94, STATUS='KEEP')
C
      RETURN
      END
