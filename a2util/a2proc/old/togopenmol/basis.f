      SUBROUTINE BASIS(IUATMS,NATOMS,ITFCT,LNP1,LNPO,NTANGM,IMEMB,
     &   NUC,NFCT,NUFCT,NANGMOM,NUMOM,ATMNAM,COORD,NPOP,NAOATM,
     &   NAOUATM,IPRINT,ISHL)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      CHARACTER*4 ATMNAM
      CHARACTER*80 XLINE
C     
      DIMENSION NUC(IUATMS),NUFCT(IUATMS),NUMOM(IUATMS),NFCT(NATOMS),
     &   NANGMOM(NATOMS),NPOP(IUATMS),NAOATM(NATOMS),NAOUATM(IUATMS),
     &   ATMNAM(IUATMS),COORD(3*NATOMS),IMEMB(NATOMS),ISHL(6)
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
      LNP1=0
      LNPO=0
      LTNP1=0
      LTNPO=0
      NTANGM=0
      DO 10 IATM=1,IUATMS
         READ(10,1110) ZNUC,IJUNK,NSHL,(ISHL(I),I=1,NSHL)
 1110    FORMAT(F20.1,8I5)
         READ(10,1115) ATMNAM(IATM),COORDX,COORDY,COORDZ
 1115    FORMAT(A4,3F20.12)
         NUFCT(IATM)=0
         NTANGMI=0
C     
         NAOTMP=0
         DO 20 I=1,NSHL
            NPT=0
            NAOT=0
            DO 21 I1=1,ISHL(I)
               READ(10,1120) NP1,NAO
               NPT=NPT+NP1
               NAOT=NAOT+NAO
               IF(I.EQ.1) NP2=1
               IF(I.EQ.2) NP2=3
               IF(I.EQ.3) NP2=6
               IF(I.EQ.4) NP2=10
               IF(I.EQ.5) NP2=15
               NAOTMP=NAOTMP+NP2*NAO
               NUFCT(IATM)=NUFCT(IATM)+NP2*NP1
               IF(I1.EQ.ISHL(I)) NTANGMI=NTANGMI+NP2
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
 10   CONTINUE
C     
      ITFCT=0
      DO 110 IATM=1,IUATMS
         DO 120 IEQATM=1,NPOP(IATM)
            ITFCT=ITFCT+NUFCT(IATM)
 120     CONTINUE
 110  CONTINUE
C     
C     Fill out NFCT, NAOATM and NMOMFCT for all atoms
      ICNT=0
      DO 1011 II=1,IUATMS
         DO 1020 IJ=1,NPOP(II)
            ICNT=ICNT+1
            NFCT(IMEMB(ICNT))=NUFCT(II)
            NANGMOM(IMEMB(ICNT))=NUMOM(II)
            NAOATM(IMEMB(ICNT))=NAOUATM(II)
 1020    CONTINUE
 1011 CONTINUE
C     
      IF(IPRINT.EQ.1)THEN
         DO 130 IATM=1,NATOMS
            WRITE(LUOUT,1170) IATM,ATMNAM(IATM),NFCT(IATM),
     &         (COORD((IATM-1)*3+J),J=1,3)
 130     CONTINUE
         WRITE(LUOUT,1130) NTANGM
         WRITE(LUOUT,1140) ITFCT
         WRITE(LUOUT,1150) LNP1
         WRITE(LUOUT,1160) LNPO
      ENDIF
 1120 FORMAT(2I5)
 1130 FORMAT(/'THE LARGEST NUMBER OF SHELLS FOR ANY ATOM IS',I3)
 1140 FORMAT('THE MOLECULE HAS',I4,' TOTAL PRIMITIVE FUNCTIONS')
 1150 FORMAT('THE LARGEST NUMBER OF PRIMITIVES IN A SINGLE SHELL IS',I3)
 1160 FORMAT('THE LARGEST NUMBER OF PRIMITIVE ORBITALS IN A SINGLE ',
     &   'SHELL IS',I3)
 1170 FORMAT(/'ATOM',I3,' IS ',A2,' WITH',I4,' PRIMITIVE FUNCTIONS'/
     &   '     ITS COORDINATES ARE ',3F15.5)
C     
      RETURN
      END
