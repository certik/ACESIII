C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      SUBROUTINE RDBAS(TITLE,IENTRY,ICONTX,LUVMOL,INTTYP,NULLST,
     &                 IATNUM,NORD,seward,cpbasis,icpun)
C
C ROUTINE READS NEW BASINF ENTRIES AND WRITES THE INFORMATION TO THE
C VMOL INPUT FILE
C
C It will scan unit 12 for the basis set and copy the definition to
C unit icpun if cpbasis=.true.
C
C THE FORMAT FOR THE ENTRIES IS AS FOLLOWS :
C
C 1. NAME OF BASIS SET.  FOLLOWS OLD BASINF CONVENTIONS, NAMELY IT
C    MUST BEGIN WITH XX:, WHERE XX IS THE ATOMIC SYMBOL.  THE NAME
C    OF THE BASIS SET BEGINS AFTER THE :.
C
C 2. THE TOTAL NUMBER OF SHELLS MAKING UP THE BASIS (NSHELL).
C
C 3. A VECTOR OF LENGTH NSHELL CONTAINING THE ANGULAR MOMENTA OF
C    THE SHELLS (LANGSH).
C
C 4. A VECTOR OF LENGTH NSHELL CONTAINING THE NUMBER OF CONTRACTED
C    FUNCTIONS IN THE SHELL (ICNTSH).
C
C 5. A VECTOR OF LENGTH NSHELL CONTAINING THE NUMBER OF PRIMITIVE
C    FUNCTIONS IN THE SHELL (IFNCSH).
C
C LOOPING OVER SHELLS
C
C1a. A VECTOR OF LENGTH IFNCSH(ISHELL) CONTAINING THE PRIMITIVE FUNCTION
C    EXPONENTS (IPRIM).
C
C1b. A MATRIX OF DIMENSION IFNCSH(ISHELL) x ICNTSH(ISHELL) CONTAINING THE
C    CONTRACTION COEFFICIENTS FOR THE SHELL.  EACH CONTRACTED FUNCTION
C    CORRESPONDS TO AN INDIVIDUAL COLUMN IN THIS MATRIX.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CJDW 9/8/97 3 lines (From Thomas Steinke).
      DOUBLE PRECISION PRT_THRESH
      CHARACTER*32 CFORMAT
      PARAMETER (PRT_THRESH = 1.0D+07 - 1.0D+00)
C      PARAMETER (MAXSHL = 20)
C      PARAMETER (NPRIMAX = 1000)
C      PARAMETER (NCNTMAX = 10)
#include "mxatms.par"
#include "baslen.par"
#include "linelen.par"
#include "baslims.par"
#include "molstrct.com"
C
      PARAMETER (NDI10 = MXATMS)
      PARAMETER (NDI9  = 250)
      parameter (lusew = 15)
      CHARACTER*(linelen) TEST
C Ajith 07/04/96
      CHARACTER*6 NULLST
      CHARACTER*(*) TITLE
      CHARACTER*(baslen) TITEL(ndi10)

      LOGICAL YESNO,seward,cpbasis
CSSS      INTEGER SHELLANG,SHELLLOC,SHELLSIZ,SHELLPRM
      DIMENSION EEXP(MXPRIM)
      DIMENSION LANGSH(MXSHEL),NINSHL(MXSHEL)
      DIMENSION IFNCSH(MXSHEL),ICNTSH(MXSHEL),COEFMT(MXPRIM*MXCORB)
C
CSSS      COMMON /MOLSTR1/ SHELLANG(100),SHELLLOC(100),SHELLSIZ(100),
CSSS     &                 SHELLPRM(100),BASISEXP(10000),BASISCNT(10000)
CSSS      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IATOM
      COMMON /TMPINF/ NSHLAT(NDI10),IATLOC(NDI10),TITEL,
     &                IOFSHL(NDI10),IOFPRM(NDI9)
      COMMON /TURBO / ITURBO,MATOM,IOFFSH
C
      INDXF(I,J,N)=I+(J-1)*N
      REWIND(12)
C
C LOCATE BASIS SET
C
1     READ(12,'(A)',END=990)TEST
      IF(INDEX(TEST,':').EQ.0)GOTO 1
      IF(TEST(1:linblnk(TEST)).NE.(TITLE(1:linblnk(TITLE))))GOTO 1
      if (cpbasis) write(icpun,'(a)') test(1:linblnk(test))
C
      READ(12,'(A)') TEST
      READ(12,*)
      READ(12,'(I3)')NSHELL
      READ(12,'((14I5))')(LANGSH(I),I=1,NSHELL)
      READ(12,'((14I5))')(ICNTSH(I),I=1,NSHELL)
      READ(12,'((14I5))')(IFNCSH(I),I=1,NSHELL)
      if (cpbasis) then
         write(icpun,'(a)') test(1:linblnk(test))
         write(icpun,*)
         write(icpun,'(I3)')NSHELL
         write(icpun,'((14I5))')(LANGSH(I),I=1,NSHELL)
         write(icpun,'((14I5))')(ICNTSH(I),I=1,NSHELL)
         write(icpun,'((14I5))')(IFNCSH(I),I=1,NSHELL)
      end if
      IDUMP=NSHELL
      IOFF=0
C
      IF(ICONTX.NE.0)THEN
C
C CHANGE CONTRACTION LEVEL IF GENBAS_X KEYWORD IS ACTIVE.  COMPLICATED
C  STUFF.
C
       ITPSHL=0
C
C DETERMINE ITPSHL (THE HIGHEST ANGULAR MOMENTUM FUNCTION IN THE BASIS)
C
       DO 10 I=1,NSHELL
        IANG=LANGSH(I)+1
        IF(ICONTX/10**(IANG-1).NE.0)ITPSHL=IANG
10     CONTINUE
C
C CALCULATE THE NUMBER OF SHELLS OF EACH ANGULAR MOMENTUM AND THE
C  NUMBER OF CONTRACTED FUNCTIONS IN EACH OF THESE SHELLS
C
       IDUMMY=ICONTX
       CALL IZERO(NINSHL,NSHELL)
       DO 11 IANGUL=1,ITPSHL
C
C NALLOW IS THE NUMBER OF CONTRACTED FUNCTIONS WHICH WILL BE USED
C FOR THIS ANGULAR MOMENTUM
C
        NALLOW=IDUMMY/10**(ITPSHL-IANGUL)
        NTOT=0
C
        DO 12 ISHELL=1,NSHELL
         IF(LANGSH(ISHELL).EQ.IANGUL-1)THEN
          NINSHL(IANGUL)=NINSHL(IANGUL)+1
          NTOT=NTOT+ICNTSH(ISHELL)
          IF(NTOT.GT.NALLOW)THEN
           NEXCESS=NTOT-NALLOW
           ICNTSH(ISHELL)=MAX(ICNTSH(ISHELL)-NEXCESS,0)
           IF(ICNTSH(ISHELL).EQ.0)THEN
            IDUMP=IDUMP-1
            NINSHL(IANGUL)=NINSHL(IANGUL)-1
           ENDIF
          ENDIF
         ENDIF
12      CONTINUE
C
        IDUMMY=IDUMMY-MIN(NTOT,NALLOW)*10**(ITPSHL-IANGUL)
C
11     CONTINUE
C

c Nevin fixed for GENBAS_X stuff, cannot use split shells and GENBAS_X
c at the same time
       NSHELL=ITPSHL

      ELSE
C
C DETERMINE ITPSHL
C
      ITPSHL=0
      DO 13 I=1,NSHELL
       ITPSHL=MAX(ITPSHL,LANGSH(I)+1)
13    CONTINUE
C
      CALL IZERO(NINSHL,ITPSHL)
      DO 14 I=1,ITPSHL
       DO 14 ISHELL=1,NSHELL
        IF(LANGSH(ISHELL).EQ.I-1) NINSHL(I)=NINSHL(I)+1
14    CONTINUE
C
      ENDIF
C
C COSMETIC STUFF.  IF HIGHEST ANGULAR MOMENTUM HAS NO SHELLS,
C DECREMENT ITPSHL BY 1.
C
      IDUMMY=0
      IOK=0
      DO 21 I=ITPSHL,1,-1
       IF(IOK.EQ.0.AND.NINSHL(I).EQ.0)THEN
        IDUMMY=IDUMMY+1
       ELSE
        IOK=1
       ENDIF
21    CONTINUE
      ITPSHL=ITPSHL-IDUMMY
C
      IF(IENTRY.EQ.0)THEN
       CALL ICOPY(NSHELL,LANGSH,1,SHELLANG(IOFFSHEL),1)
       CALL ICOPY(NSHELL,ICNTSH,1,SHELLSIZ(IOFFSHEL),1)
       CALL ICOPY(NSHELL,IFNCSH,1,SHELLPRM(IOFFSHEL),1)
       DO 313 I=1,NSHELL
        SHELLLOC(IOFFSHEL+I-1)=IATOM
313    CONTINUE
C
C STORE SHELL INFORMATION FOR THE TURBOMOLE INTERFACE
C
       IF ((ITURBO.EQ.1).AND.(INDEX(TITLE,':').NE.0)) THEN
        MATOM=MATOM+1
        IATLOC(MATOM)=MATOM
        DO 7895 I=1,MATOM
          IF (TITLE.EQ.TITEL(I)) GOTO 7896
          TITEL(MATOM)=TITLE
 7895   CONTINUE
 7896   CONTINUE
        NSHLAT(MATOM)=NSHELL
        IOFSHL(MATOM)=IOFFSHEL-1
       ENDIF

       IOFFSHEL=IOFFSHEL+NSHELL
C
        Z=FLOAT(IATNUM)
        IF(Z.GT.109.0)Z=0.0
cYAU        WRITE(LUVMOL,'(A6,F14.8,I5,(10I5))')NULLST,Z,NORD,
        WRITE(LUVMOL,'(A6,F14.8,I5,(10I5))')'      ',Z,NORD,
     &  ITPSHL,(NINSHL(I),I=1,ITPSHL)
CJDW 9/20/96. Allow for Seward.
       if(inttyp.eq.4)then
        write(lusew,'(A)')'Basis set'
        if(title(2:2).eq.':')then
          isew=1
        else
          isew=2
        endif
        write(lusew,*)title(1:isew),'..... / inline'
        Z=FLOAT(IATNUM)
        IF(Z.GT.109.0)Z=0.0
        write(lusew,'(7x,f6.1,14x,i2)')z,itpshl-1        
       ELSE IF(INTTYP.EQ.2) THEN
        WRITE(LUVMOL,'(F10.5,I5,(10I5))')Z,NORD,
     &  ITPSHL,(NINSHL(I),I=1,ITPSHL)
       ENDIF
       RETURN
      ENDIF
C
      DO 100 ISHELL=1,NSHELL
       NPRIM=IFNCSH(ISHELL)
       NCONT=ICNTSH(ISHELL)
       NJUNK=NINSHL(LANGSH(ISHELL)+1)
       READ(12,*)
       READ(12,'((5F14.6))')(EEXP(I),I=1,NPRIM)
       READ(12,*)
       if (cpbasis) then
          write(icpun,*)
          write(icpun,'((5F14.6))')(EEXP(I),I=1,NPRIM)
          write(icpun,*)
       end if
       IF(MIN(NPRIM,NCONT,NJUNK).NE.0)THEN
          WRITE(LUVMOL,'(2I5)')IFNCSH(ISHELL),ICNTSH(ISHELL)
        if(seward) then
          write(lusew,'(2i5)')nprim,ncont
        endif
       ENDIF
       IOFF0=IOFFPRIM
       DO 101 IPRIM=1,NPRIM
        LENG=MAX(1,NCONT)
        READ(12,'((7(F10.7,1X)))')(COEFMT(J),J=1,LENG)
        if (cpbasis) write(icpun,'((7(F10.7,1X)))')(COEFMT(J),J=1,LENG)
        CALL SCOPY(NCONT,COEFMT,1,BASISCNT(IOFF0),NPRIM)
        IOFF0=IOFF0+1
C
CJDW 9/8/97 5 lines (from Thomas Steinke).
        IF( EEXP(IPRIM) .GT. PRT_THRESH )THEN
          CFORMAT = '(G18.12E2,3F18.9,/,(4F18.9))'
        ELSE
          CFORMAT = '((4F18.9))'
        ENDIF
C
        IF(MIN(NJUNK,NCONT).NE.0)THEN
CJDW 9/8/97 2 lines (from Thomas Steinke).
C        WRITE(LUVMOL,'((4F18.10))')EEXP(IPRIM),(COEFMT(J),J=1,NCONT)
           WRITE(LUVMOL,CFORMAT)EEXP(IPRIM),(COEFMT(J),J=1,NCONT)
         if(seward) then
           if(iprim.eq.1)then
             do 111 isewprim=1,nprim
               write(lusew,'(f18.10)') eexp(isewprim)
 111         continue
             write(lusew,'(4f18.9)') (COEFMT(J),J=1,NCONT) 
           else
             write(lusew,'(4f18.9)') (COEFMT(J),J=1,NCONT) 
           endif
         endif
        ENDIF
101    CONTINUE
C
C STORE PRIMITIVE INFORMATION FOR THE TURBOMOLE INTERFACE
C
       IF ((ITURBO.EQ.1).AND.(INDEX(TITLE,':').NE.0)) THEN
        IOFPRM(IOFFSH+ISHELL)=IOFFPRIM-1
       ENDIF
       DO 105 IX=1,NCONT
        CALL SCOPY(NPRIM,EEXP,1,BASISEXP(IOFFPRIM),1)
        IOFFPRIM=IOFFPRIM+NPRIM
105    CONTINUE
100   CONTINUE
      IF ((ITURBO.EQ.1).AND.(INDEX(TITLE,':').NE.0))
     &   IOFFSH=IOFFSH+NSHELL
C
      if (cpbasis) write(icpun,*)
      RETURN
990   WRITE(6,5000)TITLE(1:linblnk(TITLE))
5000  FORMAT(T3,'@RDBAS-F, Basis set ',A,' not found on GENBAS.')
      CALL ERREX
C
      RETURN
      END
