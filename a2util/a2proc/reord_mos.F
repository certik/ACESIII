      SUBROUTINE REORD_MOS(IFLG,NATOMS,NBAS,NBASP,NTANGM,ICNTR,
     &                     COEFF,COEFFN,NANGMOM,NMOMAO,NAOATM,WFN)
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Logical WFN
C     
      DIMENSION COEFF(NBAS*NBASP),COEFFN(NBAS*NBASP),ICNTR(NBAS),
     $     NANGMOM(NATOMS),NMOMAO(NATOMS*NTANGM),NAOATM(NATOMS)
C     
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C     
      NBAS2=NBAS*NBASP
C
      IF (WFN .AND. IFLG.EQ.0) THEN
         CALL GETREC(20,'JOBARC','CNTERBF0',NBAS,ICNTR)
      ELSE
         IF (IFLG.EQ.0) THEN
            CALL GETREC(20,'JOBARC','CNTERBF0',NBAS,ICNTR)
            CALL GETREC(20,'JOBARC','EVECAO_A',IINTFP*NBAS2,COEFF)
         ELSE
            CALL GETREC(20,'JOBARC','EVECAO_B',IINTFP*NBAS2,COEFF)
         ENDIF
      ENDIF
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The eigenvectors in Reord_MOs"
      Call output(COEFF, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
      Write(6,*) 
      Write(6,"(a,I5)") "The value of NTANGM :", NTANGM
      Write(6,*)
      Write(6,"(a)") "The NANGMOM array"
      Write(6,"(4(1x,I3))") (NANGMOM(I), I=1, NATOMS)
      Write(6,*)
      Write(6,"(a)") "The NMOMAO array"
      Write(6,"(4(1x,I3))") (NMOMAO(I), I=1, NATOMS*NTANGM)
      Write(6,*)
#endif
C  Put higher angular momentum vectors in correct order
C
      ICNT=1
      IATMOFF=0
      DO 110 IATM=1,NATOMS
         DO 120 IANGMOM=1,NANGMOM(IATM)
            NCNT=1
            DO 125 IMOMAO=1,NMOMAO((IATM-1)*NTANGM+IANGMOM)
               IF(IANGMOM.EQ.1)THEN
                  DO 130 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=COEFF((IMO-1)*NBAS+ICNT)
 130              CONTINUE
               ENDIF
               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)
               IF(IANGMOM.EQ.2)THEN
                  DO 141 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*3+IOFF+1)
 141              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.3)THEN
                  DO 142 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*3+IOFF+2)
 142              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.4)THEN
                  DO 143 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*3+IOFF+3)
 143              CONTINUE
               ENDIF
               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)+
     $              3*NMOMAO((IATM-1)*NTANGM+2)
               IF(IANGMOM.EQ.5)THEN
                  DO 151 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+1)
 151              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.6)THEN
                  DO 152 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+2)
 152              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.7)THEN
                  DO 153 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+3)
 153              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.8)THEN
                  DO 154 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+4)
 154              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.9)THEN
                  DO 155 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+5)
 155              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.10)THEN
                  DO 156 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*6+IOFF+6)
 156              CONTINUE
               ENDIF
               IOFF=IATMOFF+NMOMAO((IATM-1)*NTANGM+1)+
     $              3*NMOMAO((IATM-1)*NTANGM+2)+
     $              6*NMOMAO((IATM-1)*NTANGM+5)
               IF(IANGMOM.EQ.11)THEN
                  DO 160 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+1)
 160              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.12)THEN
                  DO 161 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+2)
 161              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.13)THEN
                  DO 162 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+3)
 162              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.14)THEN
                  DO 163 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+4)
 163              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.15)THEN
                  DO 164 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+5)
 164              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.16)THEN
                  DO 165 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+6)
 165              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.17)THEN
                  DO 166 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+7)
 166              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.18)THEN
                  DO 167 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+8)
 167              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.19)THEN
                  DO 168 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+9)
 168              CONTINUE
               ENDIF
               IF(IANGMOM.EQ.20)THEN
                  DO 169 IMO=1,NBASP
                     COEFFN((IMO-1)*NBAS+ICNT)=
     $                    COEFF((IMO-1)*NBAS+(NCNT-1)*10+IOFF+10)
 169              CONTINUE
               ENDIF
               ICNT=ICNT+1
               NCNT=NCNT+1
 125        CONTINUE
 120     CONTINUE
         IATMOFF=IATMOFF+NAOATM(IATM)
 110  CONTINUE

#ifdef _DEBUG_LVL0        
      Write(6,*)  
      Print*, "The eigenvectors after reordering in Reord_MOs"
      Call output(COEFFN, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
      Write(6,*)  
#endif
C     
      RETURN
      END
      
      
