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

C Purpose: Produce an VMOL .MOL file based on the ZMAT file
C
C Arguments:
C     PtGrp   Symbol for the point group symmetry to use in (input)
C             making the INP file.
C     NAtms   Number of atoms (input)
C     AtNr    Atomic numbers for each center (input)
C     NUnique Number of unique internal coordinates (input)
C     geom    Cartesian coordinates (input)
c             will not be the same as Q in FINDIF
C     IStat   Error codes (output)
C             = 0  Successful
C             = 1  Fatal error from dependent
C             = 3  I/O Error on AbiFil
C             = 5  I/O Error on ZFil
C             = 7  No non-dummy centers!
C             = 8  ZFil ends after Z-matrix
C             = 10 ZFil ends after JODA control info
C     BasName Holds the basis set name given in the ZMAT file
C
C Common blocks used:
C     LUnits  Sets unit names & numbers
C     Coord
C
C Internal Variables:
C     JobTtl  Holds the job title given in the ZMAT file
C     Wrk     Scratch string to copy lines form ZMAT to ABI
C
C Dependents
C     PitPtG  Translates a schoenflies symbol into Pitzer's form
C     WrPBas  Lookup & write basis set in Pitzer style
C
C Files used:
C     ZFil   Holds the basis specification & the "rest" of INP
C     AbiFil Gets written out here (kills previous one)
C     BasFil The basis set library
C
C Limitations:
C     Only looks up bases for the geometrically unique atoms.
c
C     Extra blank lines at the end of a "truncated" ZFil will cause
C     it to act like it encountered a real error on reading the file.

#include "flags.h"

      SUBROUTINE MKVMOL(geom, PTGRP,NATMS,NUNIQUE,ZSYM,ATNR,
     $     GenBy, Scratch, ISTAT, BasName)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "io_units.par"
#include "linelen.par"
#include "baslims.par"
#include "molstrct.com"

      character*(*) BasName

      integer    lusew
      parameter (lusew  =  15)
      character*(*) SewFil
      parameter    (SewFil = 'MOLCAS.INP')

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
      Logical Opn,bad123, GENBAS
CKJW 5-24-00
      logical seward,bCpBasis,bExist,bOpened
      Integer OldLu
      Dimension Nord(2*MxAtms),Scratch(NAtms)
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
CSSS      COMMON /MOLSTR1/ SHELLANG(MXTNSH),SHELLLOC(MXTNSH),
CSSS     &                 SHELLSIZ(MXTNSH),SHELLPRM(MXTNSH),
CSSS     &                 BASISEXP(10000),BASISCNT(10000)
CSSS      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IIATOM
CSSS      COMMON /MOLSTR3/ SHELLORB(MXTNSH),SHOFFSET(MXTNSH),
CSSS     &                 PRIMORBT(MXTNSH),PROFFSET(MXTNSH)

      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/XYZIN,NWFINDIF
      COMMON /FLAGS/IFLAGS(100),IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "coord.com"
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

csb
csb      DATA BASIS /'PBS     ','STO-3G  ','DZ      ',
csb     &            'VDZ     ','3-21G   ','6-31G   ',
csb     &            '4-31G   ','TZ      ','PVTZ    ',
csb     &            '6-311G  ','DZP     ','PVDZ    ',
csb     &            'TZP     ','PVQZ    ','TZ2P    ',
csb     &            'WMR     ','6-31G*  ','6-31G** ',
csb     &            '6-311G* ','6-311G**','PV5Z    ',
csb     &            'sv      ','svp     ','dz      ',
csb     &            'dzp     ','tzp     ','tz2p    ',
csb     &            'tzpl    ','tz2pl   ','qzp     ',
csb     &            'qz2p    ','pz2p    '/

c ----------------------------------------------------------------------

C
C Initialize some things.
C

      seward = .false.

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
      CALL IGETREC(-1,'JOBARC','PASS1   ',IONE,INEWFD)
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
c
c iturbo = 1: create an interfacefile to TURBOMOLE
c
      iturbo=IFLAGS(h_IFLAGS_turbomole)
      if (iturbo.eq.1) then
        matom=0
        ioffsh=0
      endif
C
C SEE IF A GENBAS FILE IS PRESENT
C
      INQUIRE(FILE='GENBAS',EXIST=GENBAS)
      IF (.NOT.GENBAS) INQUIRE(FILE='ZMAT.BAS',EXIST=GENBAS)
C
C     Cautiously open the ZMAT file
C
      Inquire (FILE=ZFil, OPENED=Opn)
      If (.not.Opn) Open (LuZ, FILE=ZFil, STATUS='OLD')
      Rewind LuZ

C     skip the header
      bTmp = .true.
      do while (bTmp)
         read (LuZ, '(A)', ERR=8500, END=8510) JobTtl
         call parsez(JobTtl,izl)
         i = izl(1,1)
         bTmp = ((i.eq.0).or.(JobTtl(i:i).EQ.'%'))
      end do

      IF ((.NOT.XYZIN).OR.(NWFINDIF)) THEN
         i = 1
         do while (i.le.natms)
            read (LuZ,'(A)', ERR=8500, END=8510) Wrk
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
            READ (LUZ,'(A)', ERR=8500, END=8510) WRK
         END DO
      END IF
      IF (IFLAGS(h_IOPPAR_coordinates) .EQ. 0) THEN
         Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
         Do i = 1, NUnique
            Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
         End Do
      END IF
      Read (LuZ, '(A)', ERR=8500, END=8800) Wrk

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
            read(LuZ,'(A)',ERR=8500,END=8800) wrk
            bTmp = (wrk(2:2).ne.':'.and.wrk(3:3).ne.':')
         end do
c      o go back for the next read
         backspace(luz)
      end if
c YAU : end
      If (IContl.eq.0) IContl = 4
      IF (IPrnt.ge.100) THEN
         Write (LuOut, '(4X,A,I3)') 'Print Level: ', IPrnt
         Write (LuOut, '(4X,A,I3)') 'Step-taking: ', INR
         Write (LuOut, '(4X,A,I3)') 'E-Vec to follow: ',IVec
         Write (LuOut, '(4X,A,I3)') 'Die: ',IDie
         Write (LuOut, '(4X,A,I3)') 'Use curvilinear xform: ',ICurvy
         Write (LuOut, '(4X,A,I3)') 'Max Step: ',IMxStp
         Write (LuOut, '(4X,A,I3)') 'Step Contol: ', IStCrt
         Write (LuOut, '(4X,A,I3)') 'Vibrational Analysis: ', IVib
         Write (LuOut, '(4X,A,I3)') 'Recalculation of Hessian: ',IRecal
         Write (LuOut, '(4X,A,I3)') 'Convergence cutoff: ',IContl
         Write (LuOut, '(4X,A,I3)') 'Integral type: ',IntTyp
      END IF

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
               read(LuZ,'(A)',err=8500,end=8510,iostat=IOS) BasNam(i)
            else
csb 1/97 Allow arbitrary basis set names
               BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
     &                   BasName(1:linblnk(BasName))
            end if
         end if
 250  Continue
      If(IBasis.eq.0)Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
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
C     Let's find out which centers are equivalent
C
c      Call SymEqv ( NAtms, GenBy)
       do i = 1, natms
          GenBy(i) = 0
       enddo

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Count number of centers having *different* atomic numbers in this grou
C   Top half of Nord contains the number of symmetry equivalent atoms
C   for a given atomic number.
C
      Do I = 1, NAtms
         If (GenBy(i).eq.0.and.AtMass(i).gt.1.D-2) Then
            J = J + 1
            Scratch(J) = AtMass(I)
            Nord(J) = I
         End If
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

      do i = 1, 3
         ka(i) = ' '
      enddo
      NSymOp = 0

C     Now we've got (pretty much) all the information we need, so
C     we can start writing out the MOL file.
C     Begin with a very cautious opening of the file.
C
        Opn = .False. 
        Inquire (file=molfil, OPENED=Opn)
        If ( .not. Opn ) Open (LuVMol, FILE=molFil, FORM='FORMATTED',
     $     STATUS='UNKNOWN')
        Rewind LuVMol
C
C Write lines of VMOL deck.  All of this mysterious SLASK stuff, etc.
C
        idosph=1-IDFGHI
cYAU        Write(LuVMol,'(A6,4X,10I5)')Slask,IGENER,If10,If18,IECP,IECPU,
cYAU     &    IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LuVMol,'(A6,4X,10I5)')
     &     'INTGRL',IGENER,If10,If18,IECP,IECPU,
     &              IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LuVMol,'(8x,a)')
     &     ' *** ACES2 Program System (Release V0.1) *** '
        Write(LuVMol,'(A)',ERR=8300) JobTtl(:linblnk(JobTtl))
CJDW  11/2/94
        TLA = DBLE(10.0D+00**(-IFLAGS2(108)))
        Write(LuVMol,'(2I5,3A3,1X,E10.2,10X,2I5)')NonTyp,NSymop,Ka,Tla,
     &    ID3,IRstrt
        Zink=9999.0
        Zink2=3.0
        Write(LuVMol,'(2F10.2)',ERR=8300)Zink,Zink2

c   o open a basis set file for rdbas and create ZMAT.BAS if necessary
      inquire(unit=12,opened=bOpened)
      inquire(file='ZMAT.BAS',exist=bExist)
      if (bExist) then
         if (.not.bOpened) then
            open(unit=12,file='ZMAT.BAS',form='formatted',status='old')
         end if
         bCpBasis = .false.
      else
         if (.not.bOpened) then
            open(unit=12,file='GENBAS',form='formatted',status='old')
         else
            print *, '@MKVMOL: Assertion failed.'
            print *, '         ZMAT.BAS does not exist but unit 12 is ',
     &               'already open.'
            call errex
         end if
         open(unit=13,file='ZMAT.BAS',form='formatted',status='new')
         rewind(13)
         bCpBasis = .true.
      end if

C
C     Charge & multiplicity
C
c      Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
C      Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
C
C Now write the "real" input.
C
      IINFIL=0
      SHOFFSET(1)=1
      PROFFSET(1)=1
      Do 318 I = 1, Nontyp
       NStart = 1
       Do J = 1, I-1
          NStart = NStart + Nord(MxAtms+J)
       End Do
c       If(IntTyp.Eq.1)Then
c        Write(LuVMol,7000,ERR=8300)
c     &  NullSt,Float(IAtnum(Nord(NStart))),Nord(MxAtms+I)
c       Elseif(IntTyp.Eq.2)Then
c        WRITE(LUVMOL,7001,ERR=8300)
c     &  Float(IAtnum(Nord(NStart))),Nord(MxAtms+I)
c       ENDIF
c7000   Format(A6,F14.8,I5,$)
c7001   Format(F10.1,I5,$)
       SCRBAS=BASNAM(NORD(NSTART))
       IATTYP=1
       IIATOM=NORD(NSTART)
       IATOM=IATNUM(IIATOM)
       IF(IATOM.GT.2)IATTYP=2
       IF(IATOM.GT.10)IATTYP=3
       IF(GENBAS)THEN
        CALL RDBAS(SCRBAS,0,IFLAGS(64+IATTYP),LUVMOL,
     &             IntTyp,NullSt,IAtnum(Nord(NStart)),
     &             Nord(MxAtms+I),seward,.false.,13)
       ELSE
        Call WrPBas (LuBasL,BasFil,SCRBAS,
     &  LuVMol,0,1,IStat,IntTyp,Nullst,IAtnum(Nord(Nstart)),
     &  Nord(MxAtms+I))
       ENDIF

       Do J = 1, Nord(MxAtms+J)
          IBotm = -2 + 3*Nord(NStart-1+J)
          IINFIL = IINFIL + 1
          IZLOC(IINFIL) = NORD(NSTART+J-1)
          Icount = Icount + 1
          Write(LuVMol,'(A2,A1,A1,3F20.12)',ERR=8300)
     &                  ZSym(Nord(NStart)),LB,Indxx(ICOUNT),
     &                  (geom(IJ),IJ=IBotm,Ibotm+2)
       End Do
       IF(GENBAS)THEN
        CALL RDBAS(SCRBAS,1,IFLAGS(64+IATTYP),LUVMOL,
     &             ijunk,ijunk,ijunk,ijunk,seward,bCpBasis,13)
        if(inttyp.eq.4)then
          if(iflags(62).eq.0)then
            write(lusew,'(A)')'Cartesian all'
          else
            write(lusew,'(A)')'Spherical all'
          endif      
          wrk=zsym(nord(nstart))
          isew=index(wrk,' ')
          write(lusew,'(a6,3f20.12)') wrk(1:isew-1)//indxx(icount),
     &      (geom(isew),isew=ibotm,ibotm+2)
          write(lusew,'(A)')'End of Basis Set'
        endif
        SHOFFSET(I+1)=IOFFSHEL
        PROFFSET(I+1)=IOFFPRIM
        SHELLORB(I)=SHOFFSET(I+1)-SHOFFSET(I)
        PRIMORBT(I)=PROFFSET(I+1)-PROFFSET(I)
       ELSE
        Call WrPBas (LuBasL,BasFil,SCRBAS,
     &  LuVMol,1,1,IStat,ijunk,ijunk,ijunk,ijunk)
       ENDIF
C
C           Handle errors in WrPBas
C
            If (Mod(IStat,2) .eq. 1) then
               IStat = 1
               Write (LuErr, 9910)
               Close (LuZ)
               Close (LuVMol)
               Return
            ElseIf (IStat .eq. 8) then
               IStat = 3
               Write (LuErr, 9935) AbiFil, LuVMol
               Close (LuZ)
               Close (LuVMol)
               Return
            EndIf
 318    Continue
        Write(LuVMol,'(A6)',ERR=8300)'FINISH'
        if(inttyp.eq.4)then
          write(lusew,'(A)')'End of input'
        endif

c   o if we created ZMAT.BAS, then close GENBAS and connect ZMAT.BAS to unit 12
      if (bCpBasis) then
         close(12,status='keep')
         close(13,status='keep')
         open(unit=12,file='ZMAT.BAS',form='formatted',status='old')
         bCpBasis = .false.
      end if

C     Dump vector to disk which relates the atoms listed in the
C     MOL file to their position in the user supplied Z-matrix.
C
      CALL IPUTREC(20,'JOBARC','MPVMZMAT',NATOMS,IZLOC)
C
 9910 Format (' @MKDECK-F, Dependent terminated with FATAL error.')
 9935 Format (' @MKDECK-F, Dependent reports error accessing file ',
     $   A,' on unit ',I3,'.')
C
C This better not happen.
C
      If (Nontyp .eq. 0) then
         IStat = 7
         Write (LuErr, 9970)
         Close (LuZ)
         Close (LuVmol)
         Return
 9970    Format (' @MKDECK-F, No non-dummy centers in Z-matrix.')
      EndIf
C
C     The rest of this just copied from ZMAT
C
C
C     Make sure every body gets closed & stuff
C
 8000 If(IntTyp.ne.0)Close (LuVMol)
      Close (LuZ)
C
C WRITE BASIS SET INFORMATION TO JOBARC
C
      IONE=1
      NTOTPRIM=IOFFPRIM-1
      NTOTSHEL=IOFFSHEL-1
      IF(IECP.EQ.1) THEN
       CALL IPUTREC(20,'JOBARC','ATOMCHRG',NATOMS,ICRCOR)
      ELSE
       CALL IPUTREC(20,'JOBARC','ATOMCHRG',NATOMS,IATNUM)
      ENDIF
      IF(INEWFD.EQ.0)THEN
       CALL IPUTREC(20,'JOBARC','NTOTSHEL',IONE,NTOTSHEL)
       CALL IPUTREC(20,'JOBARC','NTOTPRIM',IONE,NTOTPRIM)
       CALL DPUTREC(20,'JOBARC','BASISEXP',NTOTPRIM,BASISEXP)
       CALL DPUTREC(20,'JOBARC','BASISCNT',NTOTPRIM,BASISCNT)
       CALL IPUTREC(20,'JOBARC','SHELLSIZ',NTOTSHEL,SHELLSIZ)
       CALL IPUTREC(20,'JOBARC','SHELLPRM',NTOTSHEL,SHELLPRM)
       CALL IPUTREC(20,'JOBARC','SHELLANG',NTOTSHEL,SHELLANG)
       CALL IPUTREC(20,'JOBARC','SHELLLOC',NTOTSHEL,SHELLLOC)
       CALL IPUTREC(20,'JOBARC','SHOFFSET',NONTYP,SHOFFSET)
       CALL IPUTREC(20,'JOBARC','SHELLORB',NONTYP,SHELLORB)
       CALL IPUTREC(20,'JOBARC','PROFFSET',NONTYP,PROFFSET)
       CALL IPUTREC(20,'JOBARC','PRIMORBT',NONTYP,PRIMORBT)
      ENDIF
C
      call shellinf(nontyp,natoms,ntotshel,shellorb,shoffset,shellang,
     &              shellsiz)

      Return
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Get here via an I/O error on LuAbi or LuVMol
C
 8300 IStat = 3
      If(IntTyp.ne.0)then
       Write (LuErr, 9930) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9930) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 9930 Format (' @MKDECK-F, I/O error on file ',A,'.')
C
C     Get here via an I/O error on VMLSYM
C
 9400 IStat=3
      Write (LuErr, 9930)'VMLSYM'
      Close(LuVMol)
      Close(LuZ)
      Return
C
C     Get here via an I/O error on LuZ
C
 8500 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9950) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9950) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 8510 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9955) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9955) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
C
C     If the ZFil ends after the Z-matrix or the optimization control
C     info, it probably means that JODA is being used "stand-alone"
C     and we want to allow for that possibility.
C
 8800 IStat = 8
      If(IntTyp.ne.0)then
       Write (LuOut, 9800) ZFil, MolFil
       Close (LuVMol)
      else
       Write (LuOut, 9800) ZFil, AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 8810 IStat = 10
      If(IntTyp.ne.0)then
       Write (LuOut, 9810) ZFil, MolFil
       Close (LuVMol)
      else
       Write (LuOut, 9810) ZFil, AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 9950 Format (' @MKDECK-F, I/O error on file ',A,'.')
 9955 Format (' @MKDECK-F, Premature end-of-file on ',A,'.')
 9800 Format (/1X,A,' ends after Z-matrix - cannot ',
     $   'finish making ',A,'.')
 9810 Format (/1X,A,' ends after JODA control info - cannot ',
     $   'finish making ',A,'.')
      End
