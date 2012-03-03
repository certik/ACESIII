      SUBROUTINE A2RD_BASIS_4INTGRT(IUATMS,NATOMS,ITFCT,LNP1,LNPO,
     &                              NTANGM,IMEMB,NUC,NFCT,NUFCT,
     &                              NANGMOM,NUMOM,ATMNAM,COORD,NPOP,
     &                              NAOATM,NAOUATM,NANGMOMSHL,
     &                              NCONFUNSHL,NPRIMFUNSHL,IPRINT,
     &                              ISHL,NMBERSHL,NANGMOMTSHL,
     &                              NCONFUNTSHL,NOFFSETPRM,
     &                              NOFFSETCON,NOFFSETSHL,ISHLOFF,
     &                              MAXSHELL,NUNQSHL,NTOTSHL,
     &                              IATMCHRG)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      CHARACTER*4 ATMNAM
      CHARACTER*80 XLINE
C     
      DIMENSION NUC(IUATMS),NUFCT(IUATMS),NUMOM(IUATMS),NFCT(NATOMS),
     &   NANGMOM(NATOMS),NPOP(IUATMS),NAOATM(NATOMS),NAOUATM(IUATMS),
     &   ATMNAM(IUATMS),COORD(3*NATOMS),IMEMB(NATOMS),ISHL(6), 
     &   NANGMOMSHL(IUATMS,MAXSHELL),NCONFUNSHL(IUATMS,MAXSHELL),
     &   NPRIMFUNSHL(IUATMS,MAXSHELL),NMBERSHL(NATOMS),
     &   NCONFUNTSHL(NATOMS*MAXSHELL),NANGMOMTSHL(NATOMS*MAXSHELL),
     &   NOFFSETPRM(IUATMS,MAXSHELL),NOFFSETCON(IUATMS,MAXSHELL),
     &   NOFFSETSHL(NATOMS*MAXSHELL),IATMCHRG(NATOMS),
     &   ISHLOFF(NATOMS)
C     
      COMMON /IPAR/ LUOUT
C     
C     Determine the number of functions for each atom, NFCT(IATM), the 
C     largest number of primitives in a single shell, LNP1, the largest
C     number of primitive orbitals, (number of primitive functions 
C     times the number of atomic orbitals), in a single shell, LNPO, 
C     and the total number of primitive functions, ITFCT, for the 
C     molecule.  The number of AOs for each atom NAOATM
C     
C     Open MOL for basis set information.
      OPEN(UNIT=10,FILE='MOL',FORM='FORMATTED',STATUS='OLD')
      REWIND(10)
C     

C     Read the first five lines.
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
C     
      LNP1   = 0
      LNPO   = 0
      LTNP1  = 0
      LTNPO  = 0
      NTANGM = 0
      NDUM   = 0
      IATM   = 0
      LUOUT  = 6
C
      DO IATOMS = 1, NATOMS
         IF (IATMCHRG(IATOMS) .EQ. 0) THEN
             NDUM = NDUM + 1
         Else
             IATM = IATM + 1
             IMEMB(IATM) =  IMEMB(IATM) - NDUM
         ENDIF
      ENDDO

C
      DO 10 IATM=1,IUATMS
         READ(10,1110) ZNUC,IJUNK,NSHL,(ISHL(I),I=1,NSHL)
 1110    FORMAT(F20.1,8I5)
         READ(10,1115) ATMNAM(IATM),COORDX,COORDY,COORDZ
 1115    FORMAT(A4,3F20.12)
         NUFCT(IATM)=0
         NTANGMI=0
         NUNQSHL = 0
C     
         NAOTMP=0
         DO 20 I=1,NSHL
            NPT=0
            NAOT=0
            DO 21 I1=1,ISHL(I)
               READ(10,1120) NP1,NAO
C
               NPT=NPT+NP1
               NAOT=NAOT+NAO
C
CSSS               IF(I.EQ.1) NP2=1
CSSS               IF(I.EQ.2) NP2=3
CSSS               IF(I.EQ.3) NP2=6
CSSS               IF(I.EQ.4) NP2=10
CSSS               IF(I.EQ.5) NP2=15
C
               NP2 = I*(I + 1)/2
C
               NAOTMP=NAOTMP+NP2*NAO
               NUFCT(IATM)=NUFCT(IATM)+NP2*NP1
C
C This commented if statement is wrong from the start. Ajith 01/2000
C Now this code will work for when their are multiple shells for 
C given angular momentum.
C$$$               IF(I1.EQ.ISHL(I)) NTANGMI=NTANGMI+NP2
C
                  NTANGMI=NTANGMI+NP2
                  NUNQSHL = NUNQSHL + 1 
                  NANGMOMSHL(IATM, NUNQSHL)  = NP2
                  NPRIMFUNSHL(IATM, NUNQSHL) = NP1*NP2
                  NCONFUNSHL(IATM, NUNQSHL)  = NAO*NP2
C
               NLN=(NAO-3)/4
               IF((NAO-3).GT.(NLN*4))NLN=NLN+1
               NLN=(NLN+1)*NP1
               DO 30 J=1,NLN
                  READ(10,'(A)') XLINE
 30            CONTINUE
               IF(NPT.GT.LNP1)THEN
                  IF(NPT.GT.LTNP1) LTNP1=NPT
               ENDIF
               ITMP=NPT*NAOT
               IF(ITMP.GT.LNPO)THEN
                  IF(ITMP.GT.LTNPO) LTNPO=ITMP
               ENDIF
 21         CONTINUE
 20      CONTINUE
         IF(LTNP1.GT.LNP1) LNP1=LTNP1
         IF(LTNPO.GT.LNPO) LNPO=LTNPO
         NUMOM(IATM)=NTANGMI
         IF(NTANGMI.GT.NTANGM) NTANGM=NTANGMI
         NAOUATM(IATM)=NAOTMP
         NMBERSHL(IATM) = NUNQSHL
 10   CONTINUE
C     
      ITFCT=0
      DO 110 IATM=1,IUATMS
         DO 120 IEQATM=1,NPOP(IATM)
            ITFCT=ITFCT+NUFCT(IATM)
 120     CONTINUE
 110  CONTINUE
C     
C Fill out NFCT, NAOATM and NMOMFCT for all atoms.

      ICNT=0
      DO 1011 II=1,IUATMS
C
         DO 1020 IJ = 1,NPOP(II)
C
            ICNT=ICNT+1
            NFCT(IMEMB(ICNT))=NUFCT(II)
            NANGMOM(IMEMB(ICNT))=NUMOM(II)
            NAOATM(IMEMB(ICNT))=NAOUATM(II)
            NMBERSHL(IMEMB(ICNT)) = NMBERSHL(II)
C
 1020    CONTINUE
 1011 CONTINUE
         
#ifdef _DEBUG_LVL0
      Write(6,*) "The NANGMOM array"
      Write(6, "(8i4)") (NANGMOM(I), i=1, Natoms)
      Write(6,*) "The NFCT array"
      Write(6, "(8i4)") (NFCT(I), i=1, Natoms)
      Write(6,*) "The NAOATM array"
      Write(6, "(8i4)") (NAOATM(I), i=1, Natoms)
      Write(6,*) "The NPOP array"
      Write(6, "(8i4)") (NPOP(I), i=1, Natoms)
      Write(6,*) "The NMBERSHL array"
      Write(6, "(8i4)") (NMBERSHL(I), i=1, Natoms)
#endif
C     
C Built the NCONFUNSHL, NANGMOMSHL arrays for all atoms. Also 
C get the total number of shells for all the atoms. This data
C comes handy when we try to retrive the two particle density.
C Ajith Perera 01/2001.
C

      ISHLOFF(1) = 1
      DO IATMS = 2, NATOMS
         ISHLOFF(IATMS) =  ISHLOFF(IATMS-1) + NMBERSHL(IATMS-1)
      ENDDO
#ifdef _DEBUG_LVL0      
      Write(6,*) "The offset array"
      Write(6, "(8i4)") (ISHLOFF(I), i=1, Natoms)
#endif
      NTOTSHL = 0
      ITOTATM = 0
C
      DO IATMS = 1, IUATMS
c
         DO IREDATMS = 1, NPOP(IATMS)
            ITOTATM = ITOTATM + 1
            IATMOFF = ISHLOFF(ITOTATM)

            DO IISHL = 1, NMBERSHL(IATMS)
               
               NTOTSHL = NTOTSHL + 1 

               NCONFUNTSHL(NTOTSHL) = NCONFUNSHL(IATMS, IISHL)
               NANGMOMTSHL(NTOTSHL) = NANGMOMSHL(IATMS, IISHL)
            ENDDO
         ENDDO
      ENDDO
C
      ICOUNT = 0
      DO IATMS = 1, IUATMS
C
         DO IREDATMS = 1, NPOP(IATMS)
C
            ICOUNT = ICOUNT + 1
            NMBERSHL(ICOUNT) = NMBERSHL(IATMS)
         ENDDO
      ENDDO
C
C Also bulit the offsets for primitives and contracted functions 
C on each symmetry unique atoms.
C
C$$$      CALL ICOPY(NUFCT, NUC, IUATMS)
      CALL ICOPY(IUATMS, NUFCT, 1, NUC, 1)
      NUFCT(1) = 0
      DO IATMS = 1, (IUATMS - 1)
         NUFCT(IATMS + 1) = NUC(IATMS) + NUFCT(IATMS)
      ENDDO
C
C$$$      CALL ICOPY(NAOUATM, NUC, IUATMS)
      CALL ICOPY(IUATMS, NAOUATM, 1, NUC, 1)  
      NAOUATM(1) = 0 
      DO IATMS = 1, (IUATMS - 1)
         NAOUATM(IATMS + 1) = NUC(IATMS) + NAOUATM(IATMS)
      ENDDO
C
      NOFFSETSHL(1) = 0
      DO KSHL = 1, NTOTSHL
         NOFFSETSHL(KSHL + 1) = NOFFSETSHL(KSHL) + NCONFUNTSHL(KSHL)
      ENDDO
C
C Read the atom charges back into the NUC array.
C
      CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS, NUC)
C     
      DO IATM = 1, IUATMS

         DO KSHL =1, NUNQSHL
C
            IF (IATM.EQ.1.AND.KSHL.EQ.1) THEN
C
               NOFFSETPRM(IATM,KSHL) = 0
               NOFFSETCON(IATM,KSHL) = 0
C
            ELSE IF (IATM .EQ. 1 .AND. KSHL .EQ. 2) THEN
C
               NOFFSETPRM(1, 2) = NPRIMFUNSHL(1, 1)
               NOFFSETCON(1, 2) = NCONFUNSHL(1, 1)
C
            ELSE IF (IATM .EQ. 2 .AND. KSHL .EQ. 1) THEN
C
               NOFFSETPRM(2, 1) = 0
               NOFFSETCON(2, 1) = 0 
C
            ELSE
C
               IF (KSHL .EQ. 1) THEN
                   NOFFSETPRM(IATM, KSHL) = 0
                   NOFFSETCON(IATM, KSHL) = 0
               ELSE

                   NOFFSETPRM(IATM,KSHL) = NOFFSETPRM(IATM, KSHL-1)
     &                                   + NPRIMFUNSHL(IATM, KSHL-1) 
 
                   NOFFSETCON(IATM,KSHL)= NOFFSETCON(IATM, KSHL-1)
     &                                  + NCONFUNSHL(IATM, KSHL-1) 
               ENDIF 
C
            ENDIF
         ENDDO
      ENDDO
C
#ifdef _DEBUG_LVL0
C
         WRITE(LUOUT, *) "The debug info for two particle density"
         WRITE(LUOUT, *)  "Total No, shells, Angular momentum components,
     &   No. of contracted functions per shell"
C
         DO 130 IATM=1,IUATMS
            WRITE(LUOUT,1170) IATM,ATMNAM(IATM),NFCT(IATM),
     &         (COORD((IATM-1)*3+J),J=1,3)
         IMANY_SHEL =  NMBERSHL(IATM)
         WRITE(LUOUT, *) (NANGMOMSHL(IATM,I), I= 1, IMANY_SHEL)
         WRITE(LUOUT, *) (NCONFUNSHL(IATM,I), I= 1, IMANY_SHEL)
 130     CONTINUE
         WRITE(LUOUT,1130) NTANGM
         WRITE(LUOUT,1140) ITFCT
         WRITE(LUOUT,1150) LNP1
         WRITE(LUOUT,1160) LNPO
#endif
 1120 FORMAT(2I5)
 1130 FORMAT(/'THE LARGEST NUMBER OF SHELLS FOR ANY ATOM IS',I3)
 1140 FORMAT('THE MOLECULE HAS',I4,' TOTAL PRIMITIVE FUNCTIONS')
 1150 FORMAT('THE LARGEST NUMBER OF PRIMITIVES IN A SINGLE SHELL IS',I3)
 1160 FORMAT('THE LARGEST NUMBER OF PRIMITIVE ORBITALS IN A SINGLE
     &   SHELL IS',I3)
 1170 FORMAT(/'ATOM',I3,' IS ',A2,' WITH',I4,' PRIMITIVE FUNCTIONS'/
     &   '     ITS COORDINATES ARE ',3F15.5)
C     
      RETURN
      END
