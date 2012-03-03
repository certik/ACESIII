      SUBROUTINE NODUMMY(NATMS, NTANGM, NUC, COORD, NANGMOM,
     &                   NMOMFCT,NUCTR)
C
C This routine removes dummy and ghost atoms from the molecule
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION NUC(NATMS),COORD(3*NATMS),NANGMOM(NATMS),
     &          NMOMFCT(NATMS*NTANGM),NUCTR(NATMS)
C
      ICNT=0
      DO 10 IATM=1,NATMS
         IF(NUC(IATM).NE.0)THEN
            ICNT=ICNT+1
            NUC(ICNT)=NUC(IATM)
            NUCTR(IATM)=ICNT
            NANGMOM(ICNT)=NANGMOM(IATM)
            DO 20 I=1,3
              COORD((ICNT-1)*3+I)=COORD((IATM-1)*3+I)
   20       CONTINUE
            DO 30 I=1,NTANGM
               NMOMFCT((ICNT-1)*NTANGM+I)=NMOMFCT((IATM-1)*NTANGM+I)
   30       CONTINUE
         ENDIF
   10 CONTINUE
C
      RETURN
      END
