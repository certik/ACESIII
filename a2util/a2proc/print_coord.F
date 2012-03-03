C
      SUBROUTINE PRINT_COORD(COORD, LABEL, B2ANG, NATOMS)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      CHARACTER*5 LABEL(NATOMS)
      DIMENSION COORD(NATOMS*3)
C
      Write(6,*)
      WRITE(6, 5000) 
      WRITE(6, *)
      WRITE(6, 5005)
      WRITE(6, 5010)
      WRITE(6, 5015)
      WRITE(6, 5005)
C
      Ioff =1 
      Convr2A = 1.0D0/B2Ang
C
      DO 3 I = 1, NATOMS
         WRITE(6, 5020) LABEL(I), (Convr2A*COORD(Ioff-1+K), K = 1, 3)
         Ioff = Ioff + 3
 3    CONTINUE
C
      WRITE(6, 5005)
            
 5000 FORMAT (T5, 'The normal modes corresponds to the following',
     &      /, T5,'         orientation and atomic ordering.')

 5005 FORMAT (T4, 56('-'))      
 5010 FORMAT (T20, 'Cartesian coordinates (Angs)')
 5015 FORMAT (T5, '       ', T22, 'X', T34, 'Y', T46, 'Z')
 5020 FORMAT (T5, A, T17, F10.7, T29, F10.7, T41, F10.7)
C
      RETURN
      END
