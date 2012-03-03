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
      Subroutine WrPBas (LuB,BasFil,BasNam,LuAbi,IRet,IForm,IStat,
     &                   INTTYP,NULLST,IATNUM,NORD)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: wrpbas.F,v 1.3 2010/06/30 15:58:20 ponton Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose: Read the basis library for a given basis name & write it
C          on LuAbi in the appropriate format.
C
C Arguments:
C     LuB     Unit to attach basis library to (input)
C     BasFil  Name of basis library (input)
C             NOTE: This file is opened and closed within WrPBas
C     BasNam  Name of basis set to look for (input)
C     LuAbi   Unit on which to write Pitzer-style basis info (input)
C             NOTE: This file is assumed to be open and in position
C                   on entry to WrPBas.
C     IRet    Tells if routine is to return after locating basis and
C               getting some information.
C             = 0  Just get contraction levels and highest L value.
C             = 1  Pull *and* write out basis.
C     IForm   Format to write basis set in
C             = 0  Pitzer
C             = 1  VMol
C             = 2  Cadpac
C     IStat   Completion status (output)
C             = 0  Successful
C             = 1  Named basis not found
C             = 3  Basis library file not found
C             = 5  I/O error while reading basis library
C             = 8  I/O error while writing on LuAbi
C
C Note:
C     IStat = 8 should really be a fatal error, but we don't know enough
C     to give proper error messages, so we just send it up.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer LuB, LuAbi
      Character*(*) BasFil, BasNam
C     Define the "Standard I/O" units to use throughout the code
      Parameter (LuIn  = 5)
      Parameter (LuOut = 6)
      Parameter (LuErr = 6)
C
C     For checking to open BasFil nicely
C
      Integer OldLu
      Logical IsTher, Opn
C
C     These guys are used to input & output the basis info
C
      Integer LMax, NCF, NExp
      Integer NCon(7), NPri(250), NUMSHL(250),NPRSHL(250),NINSHL(250)
      Double Precision Expnt(100), Coef(100), ARRAY(25,10)
#include "linelen.par"
      Character*(linelen) Name, Cmnt
C     Cadpac uses these & it can handle up to f fns.  Pitzer up to d.
      Character*1 LLab(7)
      double precision xyztol
      COMMON /OPTCTL/ JUNK(13),ICNTYP,JUNK2(2), xyztol
      SAVE NUMSHL,NPRSHL,NINSHL,NPRI,ISHELL,LMAX,NCF,NEXP,NAME
      Data LLab/'S','P','D','F','G','H','I'/
C
C Initialize some stuff.
C
      IsTher = .False.
      Opn = .False.
      OldLu = 0
C
C     First, carefully open the basis library
C
      Inquire (FILE=BasFil, EXIST=IsTher, OPENED=Opn, NUMBER=OldLu)
      If (.not. IsTher) then
         IStat = 3
         Write (LuErr,9930) BasFil
         Return
      EndIf
 9930 Format (' @WRPBAS-F, Basis library ',A,' not found.')
C
      If (IForm .Eq. 0 .OR. IForm .eq. 2) Then
         If (Opn) Close (OldLu)
         Open (LuB, FILE=BasFil, STATUS='OLD', FORM='FORMATTED',
     $      ERR=7000)
         Rewind LuB
      Else
         If (.not. Opn) Open (LuB, FILE=BasFil,STATUS='OLD',ERR=7000)
         If (IRet.eq.0) Rewind LuB
         If (IRet.ne.0) Then
            Backspace LuB
            NLINES=1+(LMAX+NCF)/26
            IF(MOD(LMAX+NCF,26).EQ.0)NLINES=NLINES-1
            DO 333 IBACK=1,NLINES
             Backspace LuB
333         CONTINUE
            GoTo 370
         Endif
      Endif
C
C     Now lets look...read in each basis consecutively
C
 1000 Read (LuB, '(A)', End=6000, Err=7000) Name
      Read (LuB, '(A)', End=6000, Err=7000) Cmnt
 370  Read (LuB, '(3i3)', End=6000, Err=7000) LMAX, NCF, NExp
      Read (LuB, '((26i3))', End=6000, Err=7000) (NCon(i),i=1,LMax),
     $   (NPri(i),i=1,NCF)
      If(Name(:linblnk(Name)).eq.BasNam(:linblnk(BasNam)).and.
     &   IRet.eq.0) Then
C
C Write some of this information to MOL deck and then return.
C
         ISET=0
         do i = 1, 250
            numshl(i) = 0
            nprshl(i) = 0
            ninshl(i) = 0
         enddo

         ISHELL=0
         IF(ICNTYP.EQ.0)THEN
          Do 455 LANG=1,LMAX
           NUMSHL(LANG)=NCON(LANG)
           DO 466 ISHEL=1,NUMSHL(LANG)
            ISHELL=ISHELL+1
            NINSHL(ISHELL)=NPRI(ISHELL)
            NPRSHL(ISHELL)=NPRI(ISHELL)
466        CONTINUE
455       CONTINUE
         ELSE
          Do 555 LANG=1,LMAX
           NUMSHL(LANG)=1+NCON(LANG)/8
           DO 566 ISHEL=1,NUMSHL(LANG)
            ITMP=NCON(LANG)-(ISHEL-1)*7
            ISHELL=ISHELL+1
            NINSHL(ISHELL)=MIN(7,ITMP)
            DO 567 ICOUNT=1,NINSHL(ISHELL)
             ISET=ISET+1
             NPRSHL(ISHELL)=NPRSHL(ISHELL)+NPRI(ISET)
567         CONTINUE
566        CONTINUE
555       CONTINUE
         ENDIF
CJDW 9/20/96. Allow for Seward.
CADY 5/06/04. Allow for GAMESS.
         IF(INTTYP.EQ.1.OR.INTTYP.EQ.4.OR.INTTYP.EQ.5) THEN
         WRITE(LUABI,4030)NULLST,FLOAT(IATNUM),NORD,
     &                    LMAX,(NUMSHL(I),I=1,LMAX)
         ELSE IF(INTTYP.EQ.2) THEN
         WRITE(LUABI,4031)FLOAT(IATNUM),NORD,
     &                    LMAX,(NUMSHL(I),I=1,LMAX)
         ENDIF
 4030    Format(A6,F14.8,I5,8I5)
 4031    FORMAT(F10.5,I5,8I5)
         Return
      Endif
 
      Read (LuB, '(4f18.0)', End=6000, Err=7000) (Expnt(i),i=1,NExp),
     $   (Coef(i),i=1,NExp)
C
      If ( Name(:linblnk(Name)) .ne. BasNam(:linblnk(BasNam)) )
     $   goto 1000
C
C     Now we've go to make Pitzer happy...
C
      n = 0
      nfn = 0
      IF(IForm.EQ.0)THEN
         Do 2000 l = 1, LMax
            Do 2100 ICon = 1, NCon(l)
               n = n + 1
               Write (LuAbi, 8000, Err=6100) ICon, LLab(l), NPri(n)
               Do 2200 IPri = 1, NPri(n)
                  nfn = nfn + 1
                  Write (LuAbi, 8010, Err=6100) IPri,Expnt(nfn),
     $               Coef(nfn)
 2200          Continue
 2100       Continue
 2000    Continue
 8000    Format (I5, 1X, A4, I5)
 8010    Format (I5, 1X, E14.9, E20.12)
C
C     Now we've got to make VMol happy...
C
      ElseIf (IForm .eq. 1) then
         ISTART=1
         ITHRU=0
         NROWST=1
         DO 3000 NUMSH=1,ISHELL
          CALL ZERO(ARRAY,250)
          NPRIM=NPRSHL(NUMSH)
          NSIZE=1
          IF(ICNTYP.EQ.1)NSIZE=NINSHL(NUMSH)
          WRITE(LUABI,'(2I5)')NPRIM,NSIZE
          CALL SCOPY(NPRIM,EXPNT(ISTART),1,ARRAY(1,1),1)
          IBEGIN=ISTART
          ISTART=ISTART+NPRIM
          NROWST=1
          DO 3001 I=1,NSIZE
           ITHRU=ITHRU+1
           NROW=NPRI(ITHRU)
           CALL SCOPY(NROW,COEF(IBEGIN),1,ARRAY(NROWST,I+1),1)
           IBEGIN=IBEGIN+NROW
           NROWST=NROWST+NROW
3001      CONTINUE
          DO 3002 ICOUNT=1,NPRIM
           WRITE(LUABI,3010,ERR=6100)(ARRAY(ICOUNT,J),J=1,NSIZE+1)
3002      CONTINUE
3000     CONTINUE
3010     FORMAT((4F18.9))
C
C     Cadpac is almost exactly like Pitzer, but requires fixed
C     format numbers.
C
      ElseIf (IForm .eq. 2) then
         Do 4000 l = 1, LMax
            Do 4100 ICon = 1, NCon(l)
               n = n + 1
               Write (LuAbi, 8000, Err=6100) ICon, LLab(l), NPri(n)
               Do 4200 IPri = 1, NPri(n)
                  nfn = nfn + 1
                  Write (LuAbi, 8020, Err=6100) IPri,Expnt(nfn),
     $               Coef(nfn)
 4200          Continue
 4100       Continue
 4000    Continue
 8020    Format (I5, 1X, 2F30.15)
      ENDIF
C
      Return
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     We got to the end of the file & never found a match
C
 6000 Close (LuB)
      IStat = 1
      Write (LuErr,9910) BasNam(:linblnk(BasNam)), BasFil
      Return
 9910 Format (' @WRPBAS-F, Basis set ',A,' not found in library ',A,'.')
C
C     An error occurred in writing to LuAbi
C
 6100 Close (LuB)
      IStat = 7
      Write (LuErr, 9880) LuAbi
      Return
 9880 Format (' @WRPBAS-W, I/O error writing to unit ',I3,'.')
C
C     If we get here, there's trouble - some unusual error
C
 7000 Continue
      IStat = 5
      Write (LuErr, 9950) BasFil
 9950 Format (' @WRPBAS-F, I/O error on basis library ',A,'.')
      Return
      End
