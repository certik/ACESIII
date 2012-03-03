      SUBROUTINE ADM(Q, Natoms)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (BA = 0.529177249d0)
      PARAMETER (NCOL = 6)
C
      DIMENSION Q(3*Natoms), SCR(5)
C
      IQ(I)=3*I-2

      iexit = 0

      WRITE(6,'(t3,a)') ' Interatomic distance matrix (Angstroms) '

      JBOT = 1
      NTIMES = 1 + (NATOMS-1)/5
      DO ICOUNT = 1, NTIMES
         WRITE(6,*)
         WRITE(6,144) (ICN,ICN=JBOT,MIN(NATOMS,JBOT+4))
 144     FORMAT(16X,:'[',I2,']',4(8X,:'[',I2,']'))
         DO I = JBOT, NATOMS
            ITMP = MIN(I,JBOT+4)-JBOT+1
            IF (ITMP.GT.5) THEN
               WRITE(6,*) '@ADM: Loop limit exceeded. ',I,JBOT,ITMP
            END IF
            DO J = 1, ITMP
               SCR(J) = BA * DIST(Q(IQ(I)),Q(IQ(J+JBOT-1)))
            END DO
            WRITE(6,143) I,(SCR(J),J=1,ITMP)
 143        FORMAT(T3,'[',I2,']',5(2X,F10.5))
         END DO
         JBOT = JBOT + 5
      END DO

       RETURN
       END

