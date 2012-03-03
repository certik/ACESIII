      SUBROUTINE PRIM(NATOMS,IUATMS,ITFCT,NBAS,LNP1,LNPO,NTANGM,
     &   NPOP,NFCT,NANGMOM,NMOMFCT,NMOMAO,IMEMB,NANGMOMSHL, 
     &   NPRIMFUNSHL,NCONFUNSHL,NCONFUNTSHL,NMBERSHL,ISHL2CNTR_MAP,
     &   NPRMFUNTSHL,ALPHA,PCOEFF,NAOATM,SALPHA,SPCOEF,MAXSHELL,
     &   ISHL,MAXANG)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      CHARACTER*80 XLINE
      CHARACTER*4 ATMNAM
C     
      DIMENSION NPOP(IUATMS),IMEMB(NATOMS),NFCT(NATOMS),
     &   NANGMOM(NATOMS),NMOMFCT(NATOMS*NTANGM),NAOATM(NATOMS),
     &   NMOMAO(NATOMS*NTANGM),ISHL(6),NANGMOMSHL(IUATMS,MAXSHELL),
     &   NPRIMFUNSHL(IUATMS,MAXSHELL),NCONFUNSHL(IUATMS,MAXSHELL),
     &   NCONFUNTSHL(NATOMS*MAXSHELL),NMBERSHL(IUATMS), 
     &   ISHL2CNTR_MAP(NATOMS*MAXSHELL),NPRMFUNTSHL(IUATMS*MAXSHELL)
C     
      DIMENSION PCOEFF(ITFCT*NBAS),ALPHA(ITFCT),SALPHA(LNP1),
     &   SPCOEF(LNPO)
C     
      COMMON /PAR/ PI,PIX4
      COMMON /IPAR/ LUOUT
C     
      PICNST=(0.5D+00/PI)**0.75D+00
C     
C     Read alpha and primitive coefficients from MOL file, 
C     fill in the C1 symmetry vector
      REWIND(10)
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
C     Initialize population 
      INI=1
C     
      MAXANG=0
      DO 10 IATM=1,IUATMS
         READ(10,1000) ZNUC,IJUNK,NSHL,(ISHL(I),I=1,NSHL)
 1000    FORMAT(F20.1,8I5)
         READ(10,1010) ATMNAM,COORDX,COORDY,COORDZ
 1010    FORMAT(A4,3F20.12)
C     
         IDGOFF=0
         IOFF=0
         JOFF=0
         NUNQSHL=0
         DO 28 LL=1,NSHL
            NPT=0
            NAOT=0
            DO 30 II1=1,ISHL(LL)
               READ(10,1120) NP1,NAO
               NPT=NPT+NP1
               NAOT=NAOT+NAO
               IF(LL.EQ.1) NP2=1
               IF(LL.EQ.2) NP2=3
               IF(LL.EQ.3) NP2=6
               IF(LL.EQ.4) NP2=10
               IF(LL.EQ.5) NP2=15
C
C Reinitilize the array NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL
C arrays. Ajith Perera 2001
C     
               NUNQSHL = NUNQSHL + 1 
               NANGMOMSHL(IATM, NUNQSHL)  = NP2
               NPRIMFUNSHL(IATM, NUNQSHL) = NP1
               NCONFUNSHL(IATM, NUNQSHL)  = NAO
C
               DO 31 IJ=INI,INI+NPOP(IATM)-1
                  DO 1031 IDG=1,NP2
                     NMOMFCT((IMEMB(IJ)-1)*NTANGM+IDG+IDGOFF)=NPT
                     NMOMAO((IMEMB(IJ)-1)*NTANGM+IDG+IDGOFF)=NAOT
 1031             CONTINUE
 31            CONTINUE
               IF(II1.EQ.ISHL(LL)) IDGOFF=IDGOFF+NP2
               DO 2230 I=1,LNPO
                  SPCOEF(I)=0.D+00
 2230          CONTINUE
               DO 32 I=1,NP1
                  READ(10,1060) SALPHA(I),(SPCOEF((J-1)*NP1+I),J=1,NAO)
 32            CONTINUE
C     
C     Renormalize the atomic orbitals.
C     Multiply the renormalized coefficients by the appropriate 
C     normalization constants.
C     
               DO 34 INAO=1,NAO
                  SUM=0.D+00
                  DO 36 I=1,NP1
                     DO 37 J=1,I
                        AI=SALPHA(I)
                        AJ=SALPHA(J)
                        TMP=SPCOEF((INAO-1)*NP1+I)*SPCOEF((INAO-1)*
     $                     NP1+J)*(2.0D+00*DSQRT(AI*AJ)/
     $                     (AI+AJ))**(REAL(LL)+0.5D+00)
                        SUM=SUM+TMP
                        IF(I.NE.J) SUM=SUM+TMP
 37                  CONTINUE
 36               CONTINUE
                  XNORM=1.D+00/DSQRT(SUM)
                  DO 38 I=1,NP1
                     SPCOEF((INAO-1)*NP1+I)=SPCOEF((INAO-1)*NP1+I)*
     $                  XNORM*PICNST*
     &                  (4.D+00*SALPHA(I))**(0.5D+00*REAL(LL)+0.25D+00)
 38               CONTINUE
 34            CONTINUE
C     
C     Place the alpha's and coefficients in their appropriate place in 
C     their respective matrices, ALPHA and PCOEFF.
               DO 40 IPOP=INI,INI+NPOP(IATM)-1
                  IATMOFF=0
                  JATMOFF=0
                  DO 43 ITMP=1,IMEMB(IPOP)-1
                     IATMOFF=IATMOFF+NFCT(ITMP)
                     JATMOFF=JATMOFF+NAOATM(ITMP)*ITFCT
 43               CONTINUE
                  DO 45 I=1,NP1
                     JSHOFF=0
                     DO 46 I1=1,NP2
                        ALPHA(I+IOFF+IATMOFF+(I1-1)*NP1)=SALPHA(I)
                        DO 57 J=1,NAO
                           PCOEFF(JATMOFF+(JOFF+JSHOFF)*ITFCT+I+IOFF+
     $                        IATMOFF+(I1-1)*NP1)=SPCOEF((J-1)*NP1+I)
                           JSHOFF=JSHOFF+1
 57                     CONTINUE
 46                  CONTINUE
 45               CONTINUE
 40            CONTINUE
               IOFF=IOFF+NP1*NP2
               JOFF=JOFF+NP2*NAO
 30         CONTINUE
 28      CONTINUE
         INI=INI+NPOP(IATM)
         IF(NSHL.GT.MAXANG)MAXANG=NSHL
 10   CONTINUE
C
C Reinitilize the NCONFUNTSHL array, Ajith Perer 01/2001
C
      NTOTSHL = 0
C
      DO IATMS = 1, IUATMS
C
         DO IREDATMS = 1, NPOP(IATMS)

            DO IISHL = 1, NMBERSHL(IATMS)
C            
               NTOTSHL = NTOTSHL + 1 
               NCONFUNTSHL(NTOTSHL)   = NCONFUNSHL(IATMS, IISHL)
               NPRMFUNTSHL(NTOTSHL)   = NPRIMFUNSHL(IATMS, IISHL)
               ISHL2CNTR_MAP(NTOTSHL) = IATMS
            ENDDO
         ENDDO
      ENDDO
C
 1120 FORMAT(2I5)
 1060 FORMAT(4F18.10)
      RETURN
      END
