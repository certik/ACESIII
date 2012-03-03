      SUBROUTINE A2EVAL_PRDCINT(IPRIM, JPRIM, INCRF, JNCRF, ITYPE, 
     &                          JTYPE, IPRMCOUNT, JPRMCOUNT,ICFCOUNT, 
     &                          JCFCOUNT, NTOTPRIM, NTOTCRF, MAXPRM, 
     &                          CNTMU, CNTNU, EXPS, PCOEF, CENTER, 
     &                          PRDTINT, TMP1, TMP2, TMP3)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION CNTMU(3), CNTNU(3), EXPS(NTOTPRIM), 
     &          PRDTINT(NTOTCRF, NTOTCRF), PCOEF(NTOTPRIM, NTOTCRF),
     &          CENTER(3),TMP1(MAXPRM, MAXPRM),
     &          TMP2(MAXPRM, MAXPRM),TMP3(MAXPRM, MAXPRM),
     &          AAA(27), BBB(27), DISTN(3), LMN(27), JMN(27), 
     &          CNTP(3) 
C
CFAC(9, 9)
C
      COMMON /HIGHL/ LMNVAL(3, 84), ANORM(84)
      DATA   DZERO /0.0D+00/, ONE /1.00D0/
C
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
C
C Loop over number of conttracted functions
C
      ICFINDX = ICFCOUNT
      CALL ZERO(TMP1, MAXPRM*MAXPRM)
      CALL ZERO(TMP2, MAXPRM*MAXPRM)
C
      DO ICRF = 1, INCRF
C
         ICFINDX = ICFINDX + 1 
         JCFINDX = JCFCOUNT
C 
         DO JCRF = 1, JNCRF 
C
            JCFINDX = JCFINDX + 1 
C     
C Loop over primitive functions 
C
            IPRMINDX = IPRMCOUNT
            INDX = 0
            JNDX = 0

            DO LFTPRIM = 1, IPRIM
C
               IPRMINDX = IPRMINDX + 1 
               JPRMINDX = JPRMCOUNT
               JNDX = 0
               INDX = INDX + 1 
C                     
               DO RGTPRIM = 1, JPRIM
C                  
                  JPRMINDX = JPRMINDX + 1 
                  JNDX  = JNDX + 1 
C
                  TMP1(INDX, JNDX) = 0.D0
                  TMP2(INDX, JNDX) = 0.D0
C
C We can built the product here for generate correlation hole
C
                  CALL A2BULT_PRDUCT(INDX, JNDX, ITYPE, JTYPE, 
     &                               MAXPRM, EXPS(IPRMINDX), 
     &                               EXPS(JPRMINDX), CENTER,
     &                               CNTMU, CNTNU, TMP2)
C
C Loop over primitives end here!
C
               ENDDO
            ENDDO
C
CSSS            Write(6,*) "The repulsion integral"
CSSS            CALL OUTPUT(TMP1, 1, IPRIM, 1, JPRIM, IPRIM, JPRIM, 1)
CSSS            Write(6,*) "The product integral"
CSSS            CALL OUTPUT(TMP2, 1, IPRIM, 1, JPRIM, IPRIM, JPRIM, 1)
CSSS            Write(6,*) "CONTRACTION COEFICIENTS"
CSSS            CALL OUTPUT(PCOEF, 1, NTOTPRIM, 1, NTOTCRF, NTOTPRIM,
CSSS     &                  NTOTCRF, 1)
C 
C Built the contracted functions for this shell.
C
            IOFFC = IPRMCOUNT + 1
            JOFFC = JPRMCOUNT + 1
C
C Built the contracted product functions for this shell.
C
            CALL ZERO(TMP3, MAXPRM*MAXPRM)
            CALL XGEMM('N', 'N', IPRIM, 1, JPRIM, ONE, TMP2, MAXPRM, 
     &                  PCOEF(JOFFC, JCFINDX), NTOTPRIM, DZERO, TMP3,
     &                  MAXPRM)
            CALL XGEMM('T', 'N', 1, 1, IPRIM, ONE, 
     &                  PCOEF(IOFFC, ICFINDX), NTOTPRIM, TMP3, MAXPRM,
     &                  DZERO, PRDCT, 1)
C
C$$$            Write(6,*) 
C$$$            Write(6,*) "ICFINDX, JCFINDX =", ICFINDX, JCFINDX
C$$$            Write(6,*) "REPLC, REPLD =", REPLC, REPLD

            PRDTINT(ICFINDX, JCFINDX) = PRDCT            
C
C Loop over contracted functions end here!
C
         ENDDO
      ENDDO
C
      RETURN 
      END
