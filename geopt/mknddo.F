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

C     THIS IS A MODIFICATION OF THE MKVMOL.F
C     FILE. THE PURPOSE IS TO WRITE A .MOL-LIKE
C     FILE WHICH INCLUDES ALL NON-DUMMY
C     ATOMS AND NOT JUST THE SYMMETRY UNIQUE CENTERS.
C     THIS ROUTINE IS ONLY EXECUTED WHEN THE
C     USER HAS CHOSEN AN NDDO INITIAL GUESS
C     -CARLOS TAYLOR

#include "flags.h"

      SUBROUTINE MKNDDO(geom, PTGRP,NATMS,NUNIQUE,ZSYM,ATNR,GENBY,
     $     SCRATCH, ISTAT, BasName)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "baslims.par"
#include "mxatms.par"
#include "io_units.par"
#include "linelen.par"
#include "molstrct.com"
#include "coord.com"

C SET UNIT NUMBER CARLOS TAYLOR
      parameter(lunddo=30,lunddo_input=40,ATOB=.529177d0)
      character*10 NddoFil
      character*4 CREF(3),CREFJ
      data CREF/' RHF',' UHF','ROHF'/,NddoFil/"NDDO.INPUT"/
      character*(*) BasName

      Character*(*) PtGrp, ZSym(NAtms)
      Integer AtNr(NAtms), GenBy(NAtms), IStat, AtoI,IZLoc(MxAtms)
      Character*(linelen) JobTtl, Wrk
      Character*3 Ka(3)
      Character*6 Slask,NullSt
      Character*1 Lb
      Character*1 Indxx(150)

cYAU - parse ZMAT like fetchz()
      integer izl(2,7)
      logical bTmp

csb 1/97 Hold long basis set names
#include "baslen.par"
      character*(baslen) BasNam(MxAtms),BlnkBN,ScrBas,EcpNam(MxAtms)

      dimension icrcor(mxatms)
      common /turbo / iturbo,matom,ioffsh

CJDW  11/2/94. We now use IFLAGS2(108) to set TLA.
C     Parameter (TLA = 1.D-9)
CSSS      INTEGER SHELLANG,SHELLLOC,SHELLSIZ,SHELLPRM
CSSS      INTEGER SHELLORB,SHOFFSET,PRIMORBT,PROFFSET
      double precision geom(3*mxatms)
C
      LOGICAL XYZIN,NWFINDIF
      Logical Opn,bad123
CKJW 5-24-00
      logical seward
      data seward /.false./
      Integer OldLu
      Dimension Nord(2*MxAtms),Scratch(NAtms)
C
CSSS      COMMON /MOLSTR1/ SHELLANG(100),SHELLLOC(100),SHELLSIZ(100),
CSSS     &                 SHELLPRM(100),BASISEXP(10000),BASISCNT(10000)
CSSS      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IIATOM
CSSS      COMMON /MOLSTR3/ SHELLORB(100),SHOFFSET(100),PRIMORBT(100),
CSSS     &                 PROFFSET(100)

      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/XYZIN,NWFINDIF
      COMMON /FLAGS/IFLAGS(100),IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C     Length of BasNam is fixed by the basis library
C
C
C     This length is fixed by ABINITIO's input section
C
      Data Indxx /'1','2','3','4','5','6','7','8','9','0','A'
     &,'B','C','D','E','F','G','H','I','J','K','L','M','N'
     &,'O','P','Q','R','S','T','U','V','W','X','Y','Z','a'
     &,'b','c','d','e','f','g','h','i','j','k','l','m','n'
     &,'o','p','q','r','s','t','u','v','w','x','y','z',88*'*'/

c ----------------------------------------------------------------------

C
C Initialize some things.
C
      do i=1,baslen
        BlnkBN(i:i)=' '
      end do
      ScrBas=BlnkBN

      bad123 = .false.
      Slask='INTGRL'
      NullSt='      '
      LB='#'
      J=0
      I=0
      IGENER=1
      IF10=0
      IF18=1
      ID3=0
      IONE=1
      IRSTRT=0
      IECP=IFLAGS(71)
      IECPU=0
      IPRINT=0
      INUTMP=0
      ISTATS=0
      ISTAT=0
      IOFFPRIM=1
      IOFFSHEL=1
C
      ICOUNT=0
C
C     Cautiously open the ZMAT file
C
      Inquire (FILE=ZFil, OPENED=Opn)
      If (.not.Opn) Open (LuZ, FILE=ZFil, STATUS='OLD')
      Rewind LuZ

C     skip the header
      bTmp = .true.
      do while (bTmp)
         read (LuZ, '(A)') JobTtl
         call parsez(JobTtl,izl)
         i = izl(1,1)
         bTmp = ((i.eq.0).or.(JobTtl(i:i).EQ.'%'))
      end do

      IF ((.NOT.XYZIN).OR.(NWFINDIF)) THEN
         i = 1
         do while (i.le.natms)
            read (LuZ,'(A)') Wrk
            call parsez(wrk,izl)
            if (i.eq.3) then
               izz = AtoI(wrk(izl(1,2):izl(2,2)))
               bad123 = ( izz.eq.1            .and.
     &                    zsym(1)(1:1).ne.'X' .and.
     &                    zsym(2)(1:1).ne.'X'       )
            end if
            i = i + 1
         end do
      ELSE
         DO I = 1, NATMS
            READ (LUZ,'(A)') WRK
         END DO
      END IF
      Read (LuZ,'(A)') Wrk
      IF (IFLAGS(68) .EQ. 0) THEN
         Do i = 1, NUnique
            Read (LuZ, '(A)') Wrk
         End Do
      END IF
      Read (LuZ, '(A)') Wrk

      IPrnt  = iFlags(h_IFLAGS_joda_print)
      INR    = iFlags(h_IFLAGS_opt_method)
      IContl = iFlags(h_IFLAGS_convergence)
      IVec   = iFlags(h_IFLAGS_eigenvector)
      IDie   = iFlags(h_IFLAGS_negeval)
      ICurvy = iFlags(h_IFLAGS_curvilinear)
      IStCrt = iFlags(h_IFLAGS_stp_siz_ctl)
      IMxStp = iFlags(h_IFLAGS_max_step)
      IVib   = iFlags(h_IFLAGS_vibration)
      IRecal = iFlags(h_IFLAGS_eval_hess)
      IntTyp = iFlags(h_IFLAGS_integrals)
      IDisFD = iFlags(h_IFLAGS_fd_stepsize)
      IGrdFD = iFlags(h_IFLAGS_points)
      ICNTYP = iFlags(h_IFLAGS_contraction)
      ISym   = iFlags(h_IFLAGS_symmetry)
      IBasis = iFlags(h_IFLAGS_basis)
      IDFGHI = iFlags(h_IFLAGS_spherical)
c   o if ([basis=special]) then [scroll ZMAT down to the basis set definitions]
c     WARNING: THIS IS A HACK!
c     Basis set information handling should be totally revamped.
      if (iBasis.eq.0) then
         bTmp = .true.
         do while (bTmp)
            read(LuZ,'(A)') wrk
            bTmp = (wrk(2:2).ne.':'.and.wrk(3:3).ne.':')
         end do
c      o go back for the next read
         backspace(luz)
      end if
c YAU : end
      If (IContl.eq.0) IContl = 4

C     Now read the basis sets - each line corresponds to a (non-dummy)
C     center in the Z-matrix, in the order of appearance.  We put this
C     into an array length NAtms, with blanks for dummies.  This makes
C     it easier to match the basis with the right coordinates later.
C     NOTE: blank line after basis names.
      Do 250 i = 1, NAtms
         BasNam(i) = BlnkBN
         if (AtNr(i).ne.0) then
            if (iBasis.eq.0) then
c            o read special
               read(LuZ,'(A)') BasNam(i)
            else
csb 1/97 Allow arbitrary basis set names
               BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
     &                   BasName(1:linblnk(BasName))
            end if
c        end if ([not dummy atom])
         end if
 250  Continue
      If(IBasis.eq.0)Read (LuZ,'(A)') Wrk
C ADDITION (JFS,4/90)
C Now take care of situation which occurs when 2--1--3 Z-matrix
C  specification and nonstandard basis set input is used.  In this case,
C  the basis sets for atoms #1 and #2 must be switched.
C
      If(Bad123.and.IBasis.eq.0)then
       ScrBas=BasNam(1)
       BasNam(1)=BasNam(2)
       BasNam(2)=ScrBas
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Count number of centers having *different* atomic numbers in this grou
C   Top half of Nord contains the number of symmetry equivalent atoms
C   for a given atomic number.
C
      Do I = 1, NAtms
C  MAKING THIS LOOP UNCONDITIONAL RESULTS IN ALL ATOMIC CENTERS
C  BEING WRITTEN TO FILE INSTEAD OF JUST SYMMETRY UNIQUE CENTERS
C -CARLOS TAYLOR
C         If (GenBy(i).eq.0.and.AtMass(i).gt.1.D-2) Then
            J = J + 1
            Scratch(J) = AtMass(I)
            Nord(J) = I
c         End If
      End Do
      If (IVIB.ne.2) Then
c         Call PikSr2(J,Scratch,Nord)
         NonTyp = 1
         If (J.Eq.1) Then
            Nord(MxAtms+1) = 1
         Else
            Do I = 2, J
               Dif = Dabs(Scratch(i)-Scratch(i-1))
               If (NonTyp.eq.1) IBot=0
c               If (Dif.lt.1.D-14) GoTo 314
               Nord(MxAtms+Nontyp) = I - (IBot+1)
               IBot = IBot + Nord(MxAtms+Nontyp)
               NonTyp = NonTyp + 1
 314           If (I.eq.J) Nord(MxAtms+Nontyp) = I - Ibot
            End Do
         End If
      Else If (ivib.eq.2) Then
C
C For finite difference calcs, reset this so that it writes out the
C  atoms to the MOL
C  file in exactly the same order as they appear in the ZMAT.  This
C  ordering is assumed in FINDIF!
C
         NonTyp = J
         Do I = 1, J
            Nord(MxAtms+I) = 1
         End Do
      End If

C Now get symmetry information written by SYMMETRY to file VMINF.
      If (PtGrp.EQ.'C1 ') Then
         Do I = 1, 3
            KA(I)='   '
         End Do
         NSymOp = 0
      End If
C
C     Now we've got (pretty much) all the information we need, so
C     we can start writing out the MOL file.
C     Begin with a very cautious opening of the file.
C
        Opn = .False.
        Open (LUNDDO, FILE='MOL_NDDO', FORM='FORMATTED',
     $     STATUS='UNKNOWN')
        Rewind LUNDDO
C
C Write lines of VMOL deck.  All of this mysterious SLASK stuff, etc.
C
        idosph=1-IDFGHI
cYAU        Write(LUNDDO,'(A6,4X,10I5)')Slask,IGENER,If10,If18,IECP,IECPU,
cYAU     &    IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LUNDDO,'(A6,4X,10I5)')
     &     'INTGRL',IGENER,If10,If18,IECP,IECPU,
     &              IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LUNDDO,'(8x,a)')
     &     ' *** ACES2 Program System (Release V0.1) *** '
        Write(LUNDDO,'(A)') JobTtl(:linblnk(JobTtl))
CJDW  11/2/94
        TLA = DBLE(10.0D+00**(-IFLAGS2(108)))
        Write(LUNDDO,'(2I5,3A3,1X,E10.2,10X,2I5)')NonTyp,NSymop,Ka,Tla,
     &    ID3,IRstrt
        Zink=9999.0
        Zink2=3.0
        Write(LUNDDO,'(2F10.2)')Zink,Zink2
C
C Now write the "real" input.
C
      IINFIL=0
      SHOFFSET(1)=1
      PROFFSET(1)=1
      Do 318 I = 1, Nontyp
C
C     THE FOLLOWING IF-THEN SOLVES THE PROBLEM WITH DUMMY ATOMS - CARLOS TAYLOR
C
       if(atnr(i).eq.0)goto 318
       NStart = 1
       Do J = 1, I-1
          NStart = NStart + Nord(MxAtms+J)
       End Do
       SCRBAS=BASNAM(NORD(NSTART))
       IATTYP=1
       IIATOM=NORD(NSTART)
       IATOM=IATNUM(IIATOM)
       IF(IATOM.GT.2)IATTYP=2
       IF(IATOM.GT.10)IATTYP=3
C
C     I'M NOT SURE ABOUT THE FOLLOWING CHANGE LUNDDO FOR LUVMOL
C
        CALL RDBAS(SCRBAS,0,IFLAGS(64+IATTYP),LUNDDO,
     &             IntTyp,NullSt,IAtnum(Nord(NStart)),
     &             Nord(MxAtms+I),seward,.false.,6)
       Do J = 1, Nord(MxAtms+J)
          IBotm = -2 + 3*Nord(NStart-1+J)
          IINFIL = IINFIL + 1
          IZLOC(IINFIL) = NORD(NSTART+J-1)
          Icount = Icount + 1
          Write(LUNDDO,'(A2,A1,A1,3F20.12)')
     &                  ZSym(Nord(NStart)),LB,Indxx(ICOUNT),
     &                  (geom(IJ),IJ=IBotm,Ibotm+2)
       End Do
        CALL RDBAS(SCRBAS,1,IFLAGS(64+IATTYP),LUNDDO,
     &             ijunk,ijunk,ijunk,ijunk,seward,.false.,6)
        SHOFFSET(I+1)=IOFFSHEL
        PROFFSET(I+1)=IOFFPRIM
        SHELLORB(I)=SHOFFSET(I+1)-SHOFFSET(I)
        PRIMORBT(I)=PROFFSET(I+1)-PROFFSET(I)
 318    Continue
        Write(LUNDDO,'(A6)')'FINISH'
C
C Also open the "NDDO.INPUT" file for Carlos Taylor's NDDO program. Note
C that projected NDDO is used as the initial guess for the SCF program when
C it is requested by the user (GUESS=NDDO). A. Perera, 10/2004
C
C     Determine if RHF or UHF (IREF=1 FOR RHF AND IREF=2 FOR UHF)
      IREF=IFLAGS(11)+1
      IMULT=IFLAGS(29)
      ICHG=IFLAGS(28)
      OPEN(UNIT=LUNDDO_INPUT, FILE=NddoFil,Form="Formatted",
     &     STATUS="UNKNOWN")
      WRITE(LUNDDO_INPUT,'(A40)') "'NDDO INPUT GENERATED VIA ACES2'"
      IF(IREF.EQ.1)WRITE(LUNDDO_INPUT,151)ICHG
      IF(IREF.EQ.2)WRITE(LUNDDO_INPUT,150)IMULT,ICHG
      CALL IGETREC(20, 'JOBARC', 'ATOMCHRG', NATOMS, ICRCOR)
      DO IATOMS = 1, NATOMS
         ISTART = 3*(IATOMS - 1)
         WRITE(LUNDDO_INPUT, 99) ICRCOR(IATOMS), GEOM(ISTART+ 1)*
     &                           ATOB, GEOM(ISTART+ 2)*ATOB,
     &                           GEOM(ISTART+3)*ATOB
 99      FORMAT(I2,F10.5,F10.5,F10.5)
 150     FORMAT( "'AM1 XYZ UHF MULT=",I1," CHARGE=",I1,"'")
 151     FORMAT( "'AM1 XYZ RHF CHARGE=",I1,"'")
      END DO
      CLOSE(LUNDDO)
      CLOSE(LUNDDO_INPUT)
      Close (LuZ)
C
      Return
      End
