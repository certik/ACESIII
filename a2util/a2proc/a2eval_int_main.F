      SUBROUTINE A2EVAL_INT_MAIN(MAXATMS, NATOMS, MAXSHELL, MAXPRM,
     &                           IREDUN, NSHL, NTOTPRIM, NTOTCRF, 
     &                           NOFFSETATMP, NOFFSETATMC, 
     &                           NOFFSETPRM, NOFFSETCON,NANGMOMSHL, 
     &                           NPRIMFUNSHL, NCONFUNSHL, XPOINT,
     &                           YPOINT, ZPOINT, ALPHA, PCOEF, 
     &                           COORD, PRDUTINT, TMP1, TMP2, 
     &                           TMP3)
C
C Evaluate integral of the form <Mu(r_1)Nu(r_1)> where r1 is 
C maintained as constant. Ajith Perera 01/2000 
C     
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION IP(20), FAC(9,9), IOFFST(15), COORDI(3), COORDJ(3),
     &          COORD(NATOMS*3), ALPHA(NTOTPRIM), IREDUN(MAXATMS),
     &          PCOEF(NTOTPRIM*NTOTCRF), NANGMOMSHL(MAXATMS,MAXSHELL),
     &          NPRIMFUNSHL(MAXATMS,MAXSHELL),
     &          NCONFUNSHL(MAXATMS,MAXSHELL),
     &          PRDUTINT(NTOTCRF,NTOTCRF), CENTER(3),
     &          TMP1(MAXPRM, MAXPRM), TMP2(MAXPRM, MAXPRM),
     &          TMP3(MAXPRM, MAXPRM), NSHL(MAXATMS),
     &          NOFFSETATMP(MAXATMS), NOFFSETATMC(MAXATMS), 
     &          NOFFSETPRM(MAXATMS,MAXSHELL),
     &          NOFFSETCON(MAXATMS,MAXSHELL)
C      
      CALL SETRHF(FAC, IOFFST, IP)
       
C
      IREDCNT = 0
      JREDCNT = 0
      IPRMCOUNT = 0
      JPRMCOUNT = 0
      ICFNCOUNT = 0
      JCFNCOUNT = 0
      CENTER(1) = XPOINT 
      CENTER(2) = YPOINT
      CENTER(3) = ZPOINT
C
C Loop over the symmetry unique atoms 
C      
      DO IATMS = 1, MAXATMS
C                  
         IPRMSTART =  NOFFSETATMP(IATMS)
         ICFNSTART =  NOFFSETATMC(IATMS)
         IREDATMS  =  IREDUN(IATMS)
C
C Loop over the redundent atoms 
C
         DO IRATMS = 1, IREDATMS         
C
C Loop over shells on symmetry unique atoms
C
            DO ISHL = 1, NSHL(IATMS)
C
               IREDCNT = IREDCNT + 1 
C
               DO JATMS = 1, MAXATMS
C                  
                  JPRMSTART =  NOFFSETATMP(JATMS)
                  JCFNSTART =  NOFFSETATMC(JATMS)
                  JREDATMS  =  IREDUN(JATMS)
C
                  DO JRATMS = 1, JREDATMS                 
C
                     DO JSHL = 1, NSHL(JATMS)
C
                        JREDCNT = JREDCNT + 1 
C                        
                        IANGMOM = NANGMOMSHL(IATMS, ISHL)
                        JANGMOM = NANGMOMSHL(JATMS, JSHL)
                        IPRIM   = NPRIMFUNSHL(IATMS, ISHL)
                        JPRIM   = NPRIMFUNSHL(JATMS, JSHL)
                        INCRF   = NCONFUNSHL(IATMS, ISHL)
                        JNCRF   = NCONFUNSHL(JATMS, JSHL)
                        IINTTYP = IOFFST(IANGMOM)
                        JINTTYP = IOFFST(JANGMOM)
C                        
                        ICOORD = (IATMS - 1)*3 
                        JCOORD = (JATMS - 1)*3 
C     
                        DO INDX = 1, 3
                           COORDI(INDX)  = COORD(ICOORD + INDX)
                           COORDJ(INDX)  = COORD(JCOORD + INDX)
                        ENDDO
C
C Loop over angular momentum components of each shell
C
                        IF (.NOT. ((IANGMOM*JANGMOM) .EQ. 0)) THEN
C
                           IPRMINTCNT = NOFFSETPRM(IATMS, ISHL) +
     &                                  IPRMSTART       
                           JPRMINTCNT = NOFFSETPRM(JATMS, JSHL) +
     &                                  JPRMSTART                           
                           ICFNINTCNT = NOFFSETCON(IATMS, ISHL) +
     &                                  ICFNSTART 
                           JCFNINTCNT = NOFFSETCON(JATMS, JSHL) +
     &                                  JCFNSTART 
C
                           IPRMCOUNT = IPRMINTCNT
                           ICFNCOUNT = ICFNINTCNT
                           ITYPE     = IINTTYP
C
                           DO IMOMN = 1, IANGMOM
                              ITYPE = ITYPE + 1 
                              JPRMCOUNT = JPRMINTCNT
                              JCFNCOUNT = JCFNINTCNT
                              JTYPE     = JINTTYP
C
                              DO JMOMN = 1, JANGMOM
                                 JTYPE = JTYPE + 1 
C
C Do the actual evaluation of the repulsion integral(A-Term in Eqn). Also bulit 
C the product function (B-Term in Eqn). 
C
C
                                 CALL A2EVAL_PRDCINT(IPRIM, JPRIM, 
     &                                               INCRF, JNCRF, 
     &                                               ITYPE, JTYPE, 
     &                                               IPRMCOUNT, 
     &                                               JPRMCOUNT,
     &                                               ICFNCOUNT, 
     &                                               JCFNCOUNT,
     &                                               NTOTPRIM, 
     &                                               NTOTCRF,
     &                                               MAXPRM, COORDI,
     &                                               COORDJ, ALPHA, 
     &                                               PCOEF, CENTER,
     &                                               PRDUTINT, TMP1,
     &                                               TMP2, TMP3)
C
C Angular momentum loop end here!
C     
                                 JPRMCOUNT = JPRMCOUNT + JPRIM
                                 JCFNCOUNT = JCFNCOUNT + JNCRF
                              ENDDO
                              IPRMCOUNT = IPRMCOUNT + IPRIM
                              ICFNCOUNT = ICFNCOUNT + INCRF
                           ENDDO
C     
                        ENDIF
C
C Shell loop end here!
C

                     ENDDO
                  ENDDO
               ENDDO
               JREDCNT = 0
C
            ENDDO
         ENDDO
      ENDDO
C
C$$$      V_INT = 0.0D0
CSSS      Write(6,*) "The repulsion and product integral"
C$$$      DO I =1, NTOTCRF
C$$$         DO J= 1, I
C$$$            V_INT = V_INT + REPLINT(I,J)*REPLINT(I,J)
C$$$         ENDDO
C$$$      ENDDO
CSSS      Write(6,*) "The checksum of REPLINT =", V_INT
CSSS      CALL OUTPUT(PRDUTINT, 1, NTOTCRF, 1, NTOTCRF, NTOTCRF, NTOTCRF, 1)
       
      RETURN
      END
 

