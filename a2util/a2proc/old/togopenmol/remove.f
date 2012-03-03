      SUBROUTINE REMOVE(NATOMS,NATMS,NUC,COORD,COORD2,NUC2)
C     
C     This routine removes dummy and ghost atoms from the molecule
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION NUC(NATOMS),COORD(3*NATOMS),
     $     COORD2(3*NATOMS),NUC2(NATOMS)
C     
      ICNT=0
      DO 10 IATM=1,NATOMS
         IF (NUC(IATM).NE.0) THEN
            ICNT=ICNT+1
            NUC2(ICNT)=NUC(IATM)
            DO 20 I=1,3
               COORD2((ICNT-1)*3+I)=COORD((IATM-1)*3+I)
 20         CONTINUE
         END IF
 10   CONTINUE
      NATMS=ICNT
      RETURN
      END
      
