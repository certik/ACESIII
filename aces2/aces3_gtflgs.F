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

c INPUT
c int iPrt : a flag signifying whether to dump the ASVs to stdout
c            == 0; do not print
c            != 0; print

c OUTPUT
c int iErr : an error code
c            == 0; no error
c            != 0; error <- gtflgs just dies instead of returning
c int i1-i17 : JODAFLAG record(s)
c char BasNam : the name of the basis set

#include "iachar.h"

      SUBROUTINE ACES3_GTFLGS(iPrt,iErr,
     &                  BasNam)
      IMPLICIT NONE

cJDW 3/94.
c
c Modifications have been made to extend the number of input parameters.
c /FLAGS/ has been extended to 600. In GTFLGS we have ioppar(100 + 500),
c In other subroutines, /FLAGS/ consists of IFLAGS(100) and IFLAGS2(500).
c In this way, no code in other modules which uses IFLAGS will need to be
c altered. IFLAGS2 is written to JOBARC.
c (COMMON /FLAGS2/ IFLAGS2(500)).
c
c It is intended that JFS and JG will use locations 1-100 of IFLAGS2, while
c the RJB group (and contributors other than JFS and JG) will, coordinated
c by JDW, use locations 101-200. Note that
c location N of IFLAGS2 is location 100 + N of ioppar in GTFLGS. Since
c there will be some "holes" in ioppar and nParam has been set to 600.

c Other modifications/additions :
c
c 1a. EOM-EA code (MN) :
c
c    ioppar(h_IOPPAR_ea_calc), ioppar(h_IOPPAR_ea_sym),
c    and corresponding elements of other arrays;
c    EARoot(8,2) array; EA_IRREP JOBARC record.
c    See also cMN comment lines in code.
c
c 1b. EOM-IP code (MN) :
c
c    ioppar(h_IOPPAR_ip_calc), ioppar(h_IOPPAR_ip_sym) ioppar(216) and corresponding elements of
c    other arrays; IPRoot(8,2) array; IP_IRREP JOBARC record. See also cMN
c    comment lines in code.
c
c 1c. EOMREF FLAG (MN) ioppar(h_IOPPAR_eomref) : NONE, CCSD, MBPT(2)
c     SHOULD BE SET AUTOMATICALLY BY TYPE OF CALCULATION.
c
c 2. TDHF code (HS/WJL) :
c
c    ioppar(h_IOPPAR_tdhf) used for the TDHF flag (TDHF options still currently separate
c    namelist). See also CWJL comment lines in code.
c
c 3. Extensions for number of calculation types :
c
c    NCTYPE has been increased to 50 and some of elements 26-50 have been
c    defined. Availability arrays have been extended.
c
c 4. A trap has been put in which will stop the calculation if there are no
c    *'s in the Z-matrix but METHOD has been set to something other than
c    SINGLE_POINT.
c
c 5. A trap has been put in if TDHF is attempted for other than RHF.
c
c 6. Perturbed canonical orbitals are now default for QCISD(T) gradients.
c
c 7. Numerical optimiziations and numerical hessian calculations should not
c    any longer die if analytical gradients are not available (already was
c    done for UHF/RHF numerical optimizations).
c
c 8. HFDFT code (NO) :
c
c    ioppar(204) used for the HFDFT flag which specifies which functional
c    value will be placed in JOBARC as TOTENERG
c
c 9. ABCDFULL flag added. This flag is set if value is 0 on input ('unknown)
c    It determines if abcd integrals are compressed (=2) or not (=1)
c
cJDW  7/21/94. EOM_MAXCYC flag added (element 205; element 105 of flags2)
cJDW 11/ 3/94. EOMPROP, ABCDFULL, INTGRL_TOL flags added (206-208)
cJDW  1/13/95. DAMP_TYP, DAMP_TOL, LSHF_A1, LSHF_B1 flags added (209-212)
cJDW  1/23/95. POLYRATE flag added (213).
cJDW  6/ 6/95. Add arrays iEnAva3, iGrAva3 for TWODET availability. Also
c              logical bTDCalc.
cJDW  6/ 6/95. Rename TRIM to TRMBLK.
cJDW  6/ 6/95. Add options to Ajith's EOMPROP flag; include this option
c              in logic for determining ABCDFULL.
cJDW 10/23/95. NOREORI flag added (225). This is tentative and currently
c              will only be used for POLYRATE=ON and SYMMETRY=OFF. It is
c              used in subroutine SYMMETRY.
cJDW 10/23/95. A series of modifications/additions made by MN, SRG, and AP
c              have been incorporated :
c
c              Addition of JSC_ALL option to PROPERTY ( h_IOPPAR_props)
c              Addition of P-EOMEE and BWMBPT to EXCITE ( h_IOPPAR_excite)
c              Addition of options to EA_CALC (h_IOPPAR_ea_calc) and IP_CALC
c              Modification of EOMPROP keyword options (h_IOPPAR_eomprop)
c              Addition of keyword EE_SEARCH (219)
c              Addition of keyword EOM_PRJCT (220)
c              Addition of keyword NEWVRT    (h_IOPPAR_newvrt)
c              Addition of keyword HBARABCD  (h_IOPPAR_hbarabcd)
c              Addition of keyword HBARABCI  (h_IOPPAR_hbarabci)
c              Addition of keyword KS_POT    (h_IOPPAR_ks_pot)
c
c              Addition of logic to set HBARABCD, HBARABCI
cJDW 3/26/96.  Methods CC5SD(T), CCSD-T, CC3 defined by ioppar(h_IOPPAR_calc)=31,32,33.
c              NT3EOMEE flag introduced, ioppar(224).
cJDW 7/ 9/96.  Method CCSDT-T1T2 defined by ioppar(h_IOPPAR_calc)=34.
cJDW 3/26/97.  New option GRAD_CALC (h_IOPPAR_grad_calc). Takes over role of METHOD=6, so
c              we can do transition state searches with numerical gradients.
c
cMN 7/10/97. A number of changes to incorporate MN_A3 options.
c
c              EE_SYM (h_IOPPAR_ee_sym)   : NUMBER OF ROOTS PER SYMMETRY -> 2-2-2-2/2-2-2-2
c                                                               SINGLET TRIPLET
c THIS KEYWORD CAN REPLACE 'ACES2' ESTATE_SYM -> 2/2/2/2, BOTH CAN BE USED NOW.
c
c              DIP_CALC (h_IOPPAR_dip_calc) : DOUBLE IONIZATION POTENTIALS
c              DIP_SYM (h_IOPPAR_dip_sym)  : NUMBER OF ROOTS PER SYMMETRY -> 2-2-2-2/2-2-2-2
c              DEA_CALC (h_IOPPAR_dea_calc) : DOUBLE ATTAHMENT ENERGIES
c              DEA_SYM (h_IOPPAR_dea_sym)  : NUMBER OF ROOTS PER SYMMETRY -> 2-2-2-2/2-2-2-2
c              PROGRAM (h_IOPPAR_program)  : SPECIFIES ACES2 OR MN_A3, RECOMMENDED: DEFAULT=0
c              CCR12 (233)    : CC-R12 COEFFICIENTS ARE READ IN FROM J. NOGA
c              EOMXFIELD (234): ADD FIELD IN X DIRECTION (SEE ALSO 237)
c              EOMYFIELD (235): ADD FIELD IN Y DIRECTION (SEE ALSO 237)
c              EOMZFIELD (236): ADD FIELD IN Z DIRECTION (SEE ALSO 237)
c              INSERTF (237)  : SPECIFIES POINT AT WHICH TO INSERT FIELD
c              IMEM_SIZE (239): MEMORY USED AS PSEUDO-DISK (NOT USED)
c              MAKERHF (h_IOPPAR_makerhf)  : FORCES REFERENCE TO BE RHF AS REQUIRED IN MN_A3
c              ACC_SYM (247)  : DEFINES ACTIVE SPACE IN CCSD/MBPT(2) CALCULATION
c                                  OCCUPIED/VIRTUAL, E.G. 2-1-1-0/3-2-2-2
c
c   Addition of Makerhf flag: h_IOPPAR_makerhf
c   Forces the final reference (in CCSD) to
c   be RHF. Even if it is constructed in QRHF fashion. The user should be
c   sure that at al stages in the calculation we are dealing with a
c   closed shell. Dangerous keyword!!
c
cJDW 9/16/97. Changing the way PERT_ORB keyword (option h_IOPPAR_pert_orb) works.
c             Handling case of open-shell OCC and apparent REF=RHF.
cJDW 10/30/97. Adding GLOBAL_MEM keyword for use in parallel processing
c              (but don't bet the house on it). This is option 241 (141 of
c              IFLAGS2).
cKJW 01/17/98. Added FNO_KEEP keyword option 243 (143 of IFLAGS2)
c                    FNO_POST keyword option h_IOPPAR_fno_post (144 of IFLAGS2)
c
cKJW 06/03/98. Added FNO_ACTIVE keyword option 245 (145 of IFLAGS2)
cKJW 06/09/98. Added NATURAL keyword option 246 (146 of IFLAGS2)
c
c mn 11/07/98: UNO-REF option. Three keywords
c              UNO_REF   : option h_IOPPAR_uno_ref  (logical)
c              UNO_CHARGE: option 249  ( charge of redefined state)
c              UNO_MULT  : option 250  (multiplicity of
c                          redefined state, high spin)
c AP 07/06/98: RAMAN Option.
cKJW 09/17/99: KUCHARSKI : option 252 (run Stan Kucharski's cc code)
c AP 05/12/99: SCF : Option to specify the SCF type (HF, KS ...etc)
c AP 04/2006 : Logic pertinent to vib freqency calculations diretcly 
c              followed by a geometry optimization. 

c PARAMETERS
#include "io_units.par"
#include "fnamelen.par"
#include "linelen.par"
      integer    dim_iflags,     dim_iflags2
      parameter (dim_iflags=100, dim_iflags2=500)
      INTEGER    nParam
      PARAMETER (nParam = 100 + 500)
      INTEGER    NCTYPE
      PARAMETER (NCTYPE = 50)
c     the number of ints in JODAFLAG
      INTEGER    JPARAM
      PARAMETER (JPARAM = 17)
      INTEGER    nIrMax
      PARAMETER (nIrMax = 8)
      character*4 nl_delims
      parameter  (nl_delims=',;&|')
c This can hold 8x3 64-bit integers/floats (or 192 characters)
      integer    DataSize
      parameter (DataSize = 192)

c ARGUMENTS
      integer iPrt, iErr
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17
      CHARACTER*(*) BasNam

c EXTERNAL FUNCTIONS
      INTEGER FNBLNK,LINBLNK,f_strpbrk
      CHARACTER*(linelen) TRMBLK
      CHARACTER*1 achar
      LOGICAL leq

c INTERNAL VARIABLES
      integer i, j
c     value tokens of arrays
      INTEGER nOcc(8,2),nDrop(2),iDrop(1000,2),iArr(1000)
      INTEGER EARoot(8,2), IPRoot(8,2), EERoot(8,3)
c     absolute file names from gfname
      CHARACTER*(fnamelen) FNAME
c     State Variable index returned from asv_update_kv
      integer iASV
c     logical flags
      logical bExist_JOBARC, bStar, bCCSDT3, bCompress,
     &        bTDCalc, DoABCD, DoABCI, ACES2, bMN_A3, DoQRHF, bSTEOM,
     &        EOM, bGeomOpt, bTDA, bOpenShell
      logical bOpened, bTmp, bDone, bDelLastChar, bHaveKeys, bVerbose
      logical bAutoCart, NO_AGRAD
c     character constants
      character*1 czPercent, czAsterisk, czHash, czSpace, czFirstNLChar
      character*1 czTab
c throwing implicit ints together
      integer iLength, iLine, iLines, ipt, irecal, iSpin, iUHF,
     &        nIrr, nIrrps, ndrgeo,i_havegeom, i_length, iversion
c     statistics for value tokens and line processing
      character*(linelen)  wrk, wrk2, szTmp
      character*(DataSize) szData
      integer iLastC, iOpenP, iCloseP, ValueSize, icycle

c pseudo-registers (one-time temporary integers)
      integer iTmpReg1,  iTmpReg2,  iTmpReg3,  iTmpReg4,  iTmpReg5,
     &        iTmpReg6,  iTmpReg7,  iTmpReg8,  iTmpReg9,  iTmpReg10,
     &        iTmpReg11, iTmpReg12, iTmpReg13, iTmpReg14, iTmpReg15,
     &        iTmpReg16, iTmpReg17, iTmpReg18, iTmpReg19, iTmpReg20,
     &        iTmpReg21, iTmpReg22, iTmpReg23, iTmpReg24, iTmpReg25,
     &        iTmpReg26, iTmpReg27, iTmpReg28, iTmpReg29, iTmpReg30,
     &        iTmpReg31, iTmpReg32, iTmpReg33, iTmpReg34, iTmpReg35

c COMMON BLOCKS
      INTEGER        NX,NXM6,IARCH,NCYCLE,NUNIQUE,NOPT
      COMMON /USINT/ NX,NXM6,IARCH,NCYCLE,NUNIQUE,NOPT
#define MONSTER_FLAGS /* use ioppar */
#include "jodaflags.com"
#include "machsp.com"

c ----------------------------------------------------------------------

c DECLARATIONS AND DATA STATEMENTS

c energy, gradient, and hessian availabilities for each method (technically,
c also part of some predefined structure)
      INTEGER iEnAva1(0:NCTYPE),iEnAva2(0:NCTYPE),iEnAva3(0:NCTYPE),
     &        iGrAva1(0:NCTYPE),iGrAva2(0:NCTYPE),iGrAva3(0:NCTYPE),
     &        iHsAva1(0:NCTYPE),iHsAva2(0:NCTYPE)
c RHF-UHF ENERGY CALCULATION AVAILABILITY
      DATA iEnAva1 /1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,0,0,
     &     1,1,1,1,1,
     &     1,1,0,1,0, 1,1,1,1,1, 1,1,1,1,1, 1,1,0,0,0,
     &     0,0,0,0,0/
c ROHF ENERGY CALCULATION AVAILABILITY
      DATA iEnAva2 /1, 1,1,1,1,0, 0,0,0,0,1, 0,0,1,1,0, 1,1,1,0,0,
     &     0,1,0,1,1,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 1,1,0,0,0,
     &     0,0,0,0,0 /
c TWO-DETERMINANT ENERGY CALCULATION AVAILABILITY
      DATA iEnAva3 /0, 0,0,0,0,0, 0,0,0,0,1, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 1,0,0,0,0,
     &     0,0,0,0,0 /
c RHF-UHF GRADIENT AVAILABILITY
      DATA iGrAva1 /1, 1,1,1,1,0, 0,1,1,1,1, 1,0,0,0,0, 0,0,0,0,0,
     &     1,1,1,0,0,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0/
c ROHF GRADIENT AVAILABILITY
      DATA iGrAva2 /1, 1,1,1,1,0, 0,0,0,0,1, 1,0,0,0,0, 0,0,0,0,0,
     &     0,1,0,0,0,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0/
c TWO-DETERMINANT GRADIENT CALCULATION AVAILABILITY
      DATA iGrAva3 /0, 0,0,0,0,0, 0,0,0,0,1, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 1,0,0,0,0,
     &     0,0,0,0,0 /
c RHF/UHF HESSIAN AVAILABILITY
      DATA iHsAva1 /1, 1,1,1,1,0, 0,0,1,0,1, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,1,0,0,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0/
c ROHF HESSIAN AVAILABILITY
      DATA iHsAva2 /1, 1,0,0,0,0, 0,0,0,0,1, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0,
     &     0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
     &     0,0,0,0,0/

c ----------------------------------------------------------------------

c   o define character constants
      czAsterisk = achar(_IACHAR_ASTERISK)
      czPercent  = achar(_IACHAR_PERCENT)
      czSpace    = achar(_IACHAR_SPACE)
      czHash     = achar(_IACHAR_HASH)
      czTab      = achar(_IACHAR_TAB)

c ----------------------------------------------------------------------

c   o get the external file name for JOBARC and see if it exists
      call gfname('JOBARC',fname,iLength)
      inquire(file=fname(1:iLength),exist=bExist_JOBARC)

c   o only print the logs and summaries on the first run
      call GetRec(1,'JOBARC','FIRSTRUN',1,iTmpReg1)
      if (iTmpReg1.eq.0) then
         bVerbose = .false.
         call asv_hush
      else
         bVerbose = .true.
      end if

c   o nIrrps tracks the number of irreps in all of the user's value
c     tokens. Originally, we would load COMPNIRR from JOBARC, but this
c     method fails if a low-order point group finished running in a
c     finite difference calculation.
c      if (bExist_JOBARC) then
c         call igetrec(-1,'JOBARC','COMPNIRR',1,iTmpReg1)
c      else
c         iTmpReg1 = 0
c      end if
c      nIrrps = iTmpReg1
      nIrrps = 0

c   o initialize the occupation values
      do i = 1, 16
         nOcc(i,1) = 0
      enddo

c   o initialize the dropmo values
      do i = 1, 2000
         iDrop(i,1) = 0
      enddo

      nDrop(1) = 0
      nDrop(2) = 0

c   o initialize the IP and EA irrep values
      do i = 1, 16
         EARoot(i,1) = 0
         IPRoot(i,1) = 0
      enddo

      nIrr = 8
      call iputrec(-1,'JOBARC','EA_IRREP',1,nIrr)
      call iputrec(-1,'JOBARC','IP_IRREP',1,nIrr)

      iErr = 0

c ----------------------------------------------------------------------

c LOCATE THE ACES/CRAPS NAMELIST AND COUNT THE NUMBER OF LINES IN IT

c   o open and rewind ZMAT
      inquire(file=ZFil,opened=bOpened)
      if (.not.bOpened) then
         open(unit=LUZ,file=ZFil,form='FORMATTED',status='OLD')
      end if
      rewind(luz)

c   o skip the header
      bTmp = .true.
      do while (bTmp)
         read(luz,'(a)',err=8000) wrk
         i = fnblnk(wrk)
         if (i.ne.0) then
            bTmp = (wrk(i:i).eq.czPercent).or.(wrk(i:i).eq.czHash)
         end if
      end do
c     WRK now holds the title line

c   o read the first atom line and toggle INTERNAL/CARTESIAN coordinates
      read(luz,'(a)',err=8000,end=5400) wrk
      i = 1
      do while (i.lt.80.and.
     &          (wrk(i:i).eq.czSpace.or.wrk(i:i).eq.czTab))
         i = i + 1
      end do
c     [ i=80 || wrk(i:i)!=whitespace ]
      do while (i.lt.80.and.wrk(i:i).ne.czHash.and.
     &          (wrk(i:i).ne.czSpace.and.wrk(i:i).ne.czTab))
         i = i + 1
      end do
c     [ i=80 || wrk(i:i)=whitespace || wrk(i:i)=hash ]
      do while (i.lt.80.and.
     &          (wrk(i:i).eq.czSpace.or.wrk(i:i).eq.czTab))
         i = i + 1
      end do
c     [ i=80 || wrk(i:i)!=whitespace ]
      bAutoCart = (wrk(i:i).ne.czHash.and.
     &             (wrk(i:i).ne.czSpace.and.wrk(i:i).ne.czTab))
c      if (bAutoCart) print *, '@GTFLGS: These are Cartesians.'

c   o skip the system definition but scan the lines for an asterisk
c     (which would denote a geometry optimization)
      bStar = .false.
      bTmp  = .true.
      do while (bTmp)
         read(luz,'(a)',err=8000,end=5400) wrk
         i = fnblnk(wrk)
         if (i.ne.0) then
            if (wrk(i:i).eq.czAsterisk) then
               bTmp = .false.
            else
c            o stop at the first * or # (thus ignoring *'s in comments)
               i = f_strpbrk(wrk,czAsterisk//czHash)
               if (i.ne.0) then
                  bStar = bStar .or. (wrk(i:i).eq.czAsterisk)
               end if
            end if
         end if
      end do
c     WRK now holds line 1 of the first namelist

c   o skip to the ACES2/CRAPS namelist
      bTmp = .true.
      do while (bTmp)
         if ( (wrk(i:(i+5)).eq.'*ACES2') .or.
     &        (wrk(i:(i+5)).eq.'*CRAPS')      ) then
            bTmp = .false.
         else
            i = 0
            do while (i.eq.0)
               read(luz,'(a)',err=8000,end=5400) wrk
               i = fnblnk(wrk)
               if (i.ne.0) then
                  if (wrk(i:i).ne.czAsterisk) i = 0
               end if
            end do
         end if
      end do

c   o While we're counting lines, record the first and last characters
c     for future reconditioning.
      czFirstNLChar = czSpace
      bDelLastChar  = .false.

c   o count the number of lines with keywords on them
c     There are three ways to terminate the ACES2/CRAPS namelist:
c      - there is an unpaired close parenthesis, ')'
c      - an asterisk is the first non-blank character on the line
c      - the end-of-file is reached
      iLines = 1
      bDone  = .false.
      do while (.not.bDone)

c      o recondition the test string and point to the first char
c        (both in the namelist, with czFirstNLChar, and in the line, with i)
         iLastC = index(wrk,czHash)
         if (iLastC.eq.0) then
            wrk2 = trmblk(wrk)
            wrk  = wrk2
            i    = fnblnk(wrk)
         else
c         o pure comments will not terminate the namelist
            i = 1
            if (iLastC.eq.1) then
               wrk = ' '
            else
               wrk2 = trmblk(wrk(1:iLastC-1))
               wrk  = wrk2
            end if
         end if
         if (czFirstNLChar.eq.czSpace) then
            if (iLines.eq.1) then
               j = 6+fnblnk(wrk(7:))
               if (j.ne.6) czFirstNLChar = wrk(j:j)
            else
               j = fnblnk(wrk)
               if (j.ne.0) czFirstNLChar = wrk(j:j)
            end if
         end if

c      o test for a (totally) blank line or a new namelist after the first line
c        (EOF automatically jumps)
         if (i.ne.0) then
            bTmp = (wrk(i:i).eq.czAsterisk)
         else
            bTmp = .false.
         end if
         if ((i.eq.0).or.((iLines.ne.1).and.bTmp)) then
c         o read too far -> jump back and quit counting
            backspace(luz)
            iLines = iLines - 1
            bDone  = .true.
c         o Did we accidentally record an asterisk for the first char of
c           an empty namelist?
            if (czFirstNLChar.eq.czAsterisk) czFirstNLChar = czSpace

         else
c         o wrk has no comments and may be empty (denoting a pure comment
c           in the namelist, which won't terminate parsing)
            iCloseP = index(wrk,')')
            if (iCloseP.eq.0) then
c            o no ')' -> keep going
c              (NOTE: we will check for multiple, unmatched '(' during the
c                     key-value tokenization)
               read(luz,'(a)',end=8012) wrk
               iLines = iLines + 1

c           else if (iCloseP.ne.0) then
            else
c            o measure the parenthetical level (j) from the beginning (i)
               i = 0
               if (iLines.eq.1) then
                  j = 0
               else
c               o This may seem strange, but we have to convince the logic
c                 after the parenthetical count that we are still in the
c                 namelist. If we encounter a closing ')', then we will test
c                 czFirstNLChar to see if it is '('.
                  j = 1
               end if

               iTmpReg2 = f_strpbrk(wrk,'()')
               do while (iTmpReg2.ne.0)

                  i = i + iTmpReg2
                  if (wrk(i:i).eq.'(') then
                     j = j + 1
                  else
                     j = j - 1
                  end if

                  if (i.lt.len(wrk)) then
                     iTmpReg2 = f_strpbrk(wrk(i+1:),'()')
                  else
                     iTmpReg2 = 0
                  end if

c              end do while ([have more parentheses])
               end do

c            o no more parentheses
               if (j.eq.0) then
c               o we should be done with the namelist, but is there any data
c                 after the last ')'?
                  if (i.lt.len(wrk)) then
                     iTmpReg3 = fnblnk(wrk(i+1:))
                     if (iTmpReg3.ne.0) then
                        print *, '@GTFLGS: ERROR - ',
     &         'data found after the terminal ")" in the ACES namelist:'
                        print *, wrk
                        call errex
                     end if
                  end if
c               o since a final ')' was found, make sure the namelist was
c                 started with a '('
                  if (czFirstNLChar.ne.'(') then
                     print *, '@GTFLGS: ERROR - ")" terminates the ',
     &                        'namelist but "(" does not initiate it.'
                     call errex
                  end if
                  bDelLastChar = .true.
                  bDone = .true.

c              else if (j.ne.0) then
               else

                  if (j.eq.1) then
c                  o we should still be somewhere in the namelist
                     read(luz,'(a)',end=8012) wrk
                     iLines = iLines + 1
                  else
                     print *, '@GTFLGS: Error - ',
     &                        'open parentheses in the ACES namelist:'
                     print *, wrk
                     call errex
                  end if

c              end if (j.eq.0)
               end if

c           end if (iCloseP.eq.0)
            end if

c        end if ([blank line].or.[new namelist])
         end if
c     end do while (.not.bDone)
      end do

c   o EOF end jump
 8012 continue

c   o If the namelist is empty, then stop.
      if ((iLines.eq.1).and.(czFirstNLChar.eq.czSpace)) then
         print *, '@GTFLGS: ACES will not run with an empty namelist.'
         call errex
      end if

c ----------------------------------------------------------------------

c PARSE THE ACES NAMELIST

c   o return to the beginning of the namelist
c      print *, 'DEBUG: Lines in namelist = ',iLines
      do i = 1, iLines
         backspace(luz)
      end do

c   o print the ASV registration header
      if (bVerbose) then
         print '(/)'
         print *, '    ',
     &   '                 ACES STATE VARIABLE REGISTRATION LOG'
         print *, '    ', ('-',iTmpReg4=1,70)
      end if

c   o initialize flags to default values
c      call asv_update_kv(achar(0),iASV,wrk,0)
      call init_flags()
      bHaveKeys = .false.

c   o loop over the number of lines
      do iLine = 1, iLines
         read(luz,'(a)') wrk

c      o recondition the token string
         iLastC = index(wrk,czHash)
         if (iLastC.eq.1) then
            wrk = ' '
         else
            if (iLastC.eq.0) then
               wrk2 = trmblk(wrk)
            else
               wrk2 = trmblk(wrk(1:iLastC-1))
            end if
            if (iLine.eq.1) then
c            o throw out the name of the namelist
               wrk  = wrk2(7:)
               wrk2 = trmblk(wrk)
            end if
c         o throw out an opening '('
c           The only time an open parenthesis should appear in the namelist
c           before a key token is in the first line, so this condition should
c           be valid.
            if (f_strpbrk('(',wrk2(1:1)).ne.0) wrk2(1:1) = czSpace
            if ((iLine.eq.iLines).and.bDelLastChar) then
               iTmpReg5 = linblnk(wrk2)
               wrk2(iTmpReg5:iTmpReg5) = czSpace
            end if
            wrk = trmblk(wrk2)
         end if
         wrk2 = wrk

c      print *, 'DEBUG: Attempting to parse'
c      print *, wrk

c      o loop over boolean and key-value tokens
c        (ipt is the roving pointer. It will never point to a delimiter.)
         ipt = fnblnk(wrk)
         bDone = (ipt.eq.0)
         do while (.not.bDone)

            if (ipt.gt.len(wrk)) then
               bDone = .true.
               goto 101
            end if

c         o get off a delimiter and skip empty sets (like ",, ,  ,")
            iTmpReg6 = f_strpbrk(nl_delims,wrk(ipt:ipt))
            do while (iTmpReg6.ne.0)
               ipt = ipt + 1
               if (ipt.gt.len(wrk)) then
                  bDone = .true.
                  goto 101
               end if
               iTmpReg6 = fnblnk(wrk(ipt:))
               if (iTmpReg6.ne.0) then
                  ipt  = ipt-1 + iTmpReg6
                  iTmpReg6 = f_strpbrk(nl_delims,wrk(ipt:ipt))
               end if
            end do

c         o find the next delimiter (boolean) or '=' (kvpair)
            iTmpReg7 = f_strpbrk(wrk(ipt:),nl_delims//'=')
            if (iTmpReg7.eq.0) then
c            o check for a final boolean key token
               iTmpReg8 = fnblnk(wrk(ipt:))
               if (iTmpReg8.ne.0) then
                  iLastC = linblnk(wrk)
                  call asv_update_kv(wrk(ipt:iLastC)//achar(0),iASV,
     &                               szData,0)
                  bHaveKeys = .true.
               end if
               bDone = .true.
               goto 101
            else
               iTmpReg9 = ipt-1 + iTmpReg7
               if (wrk(iTmpReg9:iTmpReg9).eq.'=') then
c               o find the end of the value token (beware empty keys)
                  if (iTmpReg9.eq.ipt) then
                     print *, '@GTFLGS: Missing key token.'
                     print *, wrk
                     call errex
                  end if
                  if (iTmpReg9.eq.len(wrk)) then
                     print *, '@GTFLGS: Missing value token.'
                     print *, wrk
                     call errex
                  end if
                  j = f_strpbrk(wrk(iTmpReg9+1:),nl_delims)
                  if (j.eq.0) then
                     iLastC = linblnk(wrk)
                  else
                     iLastC = iTmpReg9+j-1
                  end if
c               o pass only the key-value tokens to the ASV updater
                  szTmp = wrk(ipt:iLastC)//achar(0)
                  ValueSize = DataSize
cOLD                  call asv_update_kv(szTmp,iASV,szData,ValueSize)
                  call asv_update_kv(wrk(ipt:iLastC)//achar(0),iASV,
     &                               szData,ValueSize)
c               o point to the next key token
                  ipt = iLastC + 2
               else
c               o process the boolean and point to the next key token
                  call asv_update_kv(wrk(ipt:iTmpReg9-1)//achar(0),iASV,
     &                               szData,0)
                  bHaveKeys = .true.
                  ipt = iTmpReg9 + 1
                  goto 101
               end if
            end if

c         o What do we do if the ASV is not recognized: die or skip?
            if (iASV.eq.0) call errex
            bHaveKeys = .true.
c            if (iASV.eq.0) goto 101

c         o Now, process the ASV if it is of type string, array, or double
c           This is not at all efficient in the sense that we if/endif
c           for EVERY ASV, but since Fortran has no switch/case construct,
c           we will deal with self-contained if/endif blocks. You could
c           drop a 'goto 101' at the end of each block, but this might
c           impair performance (as if it counts at this point).

c         o Note: A few ASVs take data up to the maximum number of allowable
c           irreps. OCCUPATION and IP_SYM are two such ASVs. nIrrps is
c           initialized to zero. The first irrep-defining key that is read
c           will set this value. All others will test against it and pass
c           or fail.

c         o BASIS
            if (iASV.eq.h_IOPPAR_basis) then
               if (leq(szData(1:ValueSize),'SPECIAL')) then
                  ioppar(h_IOPPAR_basis) = 0
               else
                  BasNam = szData(1:ValueSize)
                  ioppar(h_IOPPAR_basis) = 1
                  iTmpReg10 = 1+ValueSize
                  call iputrec(1,'JOBARC','BASNAMLN',1,iTmpReg10)
                  call putcrec(1,'JOBARC','BASISNAM',iTmpReg10,BasNam)
               end if
               goto 101
            end if

c YAU - add a conditional for SUBGROUP that sets nIrrps

c         o OCCUPATION
            if (iASV.eq.h_IOPPAR_occ) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             nOcc,8,2)
               ioppar(h_IOPPAR_occ) = iSpin
               goto 101
            end if

c         o MN_A3 symmetry arrays
            if (iASV.eq.h_IOPPAR_ea_sym) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             EARoot,8,2)
               call PutRec(-1,'JOBARC','EA_IRREP',1,nIrrps)
               ioppar(h_IOPPAR_ea_sym) = iSpin
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_ip_sym) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             IPRoot,8,2)
               call PutRec(-1,'JOBARC','IP_IRREP',1,nIrrps)
               ioppar(h_IOPPAR_ip_sym) = iSpin
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_ee_sym) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             EERoot,8,3)
               call iputrec(20,'JOBARC','EESYM_A ',8,EERoot(1,1))
               call iputrec(20,'JOBARC','EESYM_B ',8,EERoot(1,2))
               call iputrec(20,'JOBARC','EESYM_C ',8,EERoot(1,3))
               call iputrec(20,'JOBARC','ICOUNT  ',1,nIrrps)
               call iputrec(20,'JOBARC','EESYMINF',nIrrps,EERoot(1,1))
               ioppar(h_IOPPAR_ee_sym) = iSpin
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_dip_sym) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             EERoot,8,3)
               call iputrec(20,'JOBARC','DIPSYMA ',8,EERoot(1,1))
               call iputrec(20,'JOBARC','DIPSYMB ',8,EERoot(1,2))
               call iputrec(20,'JOBARC','DIPSYMC ',8,EERoot(1,3))
               ioppar(h_IOPPAR_dip_sym) = iSpin
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_dea_sym) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             EERoot,8,3)
               call iputrec(20,'JOBARC','DEASYMA ',8,EERoot(1,1))
               call iputrec(20,'JOBARC','DEASYMB ',8,EERoot(1,2))
               call iputrec(20,'JOBARC','DEASYMC ',8,EERoot(1,3))
               ioppar(h_IOPPAR_dea_sym) = iSpin
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_acc_sym) then
               call parse_irp_iarr(szData(1:ValueSize),nIrrps,iSpin,
     &                             EERoot,8,2)
               call iputrec(20,'JOBARC','CCSYM_O ',8,EERoot(1,1))
               call iputrec(20,'JOBARC','CCSYM_V ',8,EERoot(1,2))
               ioppar(h_IOPPAR_acc_sym) = iSpin
               goto 101
            end if

c         o DROPMO
            if (iASV.eq.h_IOPPAR_dropmo) then
               call parse_set_iarr(szData(1:ValueSize),j,iDrop,2000)
               nDrop(1) = j
               ioppar(h_IOPPAR_dropmo) = j
               goto 101
            end if

c         o FD_IRREPS
            if (iASV.eq.h_IOPPAR_fd_irreps) then
               call parse_set_iarr(szData(1:ValueSize),j,iArr,1000)
               call iputrec(20,'JOBARC','NFDIRREP',1,j)
               call iputrec(20,'JOBARC','FDIRREP ',j,iArr)
               ioppar(h_IOPPAR_fd_irreps) = j
               goto 101
            end if

c         o ESTATE_SYM
            if (iASV.eq.h_IOPPAR_estate_sym) then
               call parse_set_iarr(szData(1:ValueSize),j,iArr,8)
               call iputrec(20,'JOBARC','ICOUNT  ',1,j)
               call iputrec(20,'JOBARC','EESYMINF',j,iArr)
               if (j.lt.8) then
                  do i = 1, 8-j
                     iArr(j+i) = 0
                  enddo
               endif
               call iputrec(20,'JOBARC','EESYM_A ',8,iArr)
               do i = 1, 8
                  iArr(i) = 0
               enddo
               call iputrec(20,'JOBARC','EESYM_B ',8,iArr)
               call iputrec(20,'JOBARC','EESYM_C ',8,iArr)
               ioppar(h_IOPPAR_estate_sym) = j
               goto 101
            end if

c         o QRHF orbital arrays
            if (iASV.eq.h_IOPPAR_qrhf_gen) then
               call parse_set_iarr(szData(1:ValueSize),j,iArr,1000)
               call iputrec(20,'JOBARC','QRHFTOT ',1,j)
               call iputrec(20,'JOBARC','QRHFIRR ',j,iArr)
               ioppar(h_IOPPAR_qrhf_gen) = j
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_qrhf_orb) then
               call parse_set_iarr(szData(1:ValueSize),j,iArr,1000)
               call iputrec(20,'JOBARC','QRHFLOC ',j,iArr)
               ioppar(h_IOPPAR_qrhf_orb) = j
               goto 101
            end if
            if (iASV.eq.h_IOPPAR_qrhf_spin) then
               call parse_set_iarr(szData(1:ValueSize),j,iArr,1000)
               call iputrec(20,'JOBARC','QRHFSPN ',j,iArr)
               ioppar(h_IOPPAR_qrhf_spin) = j
               goto 101
            end if

  101       continue
c        end do while (.not.bDone) <- more bools or kvpairs left on line
         end do

c     end do iLine = 1, iLines
      end do

c   o print the ASV registration footer
      if (bVerbose) print *, '    ', ('-',iTmpReg10=1,70)

c   o If the namelist is empty, then stop.
      if (.not.bHaveKeys) then
         print *, '@GTFLGS: ACES will not run with an empty namelist.'
         call errex
      end if

c ----------------------------------------------------------------------

c IMMEDIATE DEATH!!!

c   o If the reference is not RHF, make sure OCC has both spins defined.
c     (Would it be so unreasonable to simply copy the alpha occs into beta?)
      if (ioppar(h_IOPPAR_ref).ne.0.and.
     &    ioppar(h_IOPPAR_occ).eq.1     ) then
          print *, '@GTFLGS: OCCUPATION must define both alpha and beta'
          print *, '         spins, since the reference is not spin ',
     &             'restricted.'
          close(unit=luz,status='KEEP')
          call errex
      end if

c   o two QRHF keys are dependent on qrhf_gen
      if (ioppar(h_IOPPAR_qrhf_orb).ne.0                        .and.
     &    ioppar(h_IOPPAR_qrhf_orb).ne.ioppar(h_IOPPAR_qrhf_gen)
     &   ) then
         print *, '@GTFLGS: The number of orbitals listed in QRHF_ORB',
     &            ' is not the same as'
         print *, '         the number of orbitals listed in QRHF_GEN.'
         call errex
      end if
      if (ioppar(h_IOPPAR_qrhf_spin).ne.0                        .and.
     &    ioppar(h_IOPPAR_qrhf_spin).ne.ioppar(h_IOPPAR_qrhf_gen)
     &   ) then
         print *, '@GTFLGS: The number of orbitals listed in QRHF_SPIN',
     &            ' is not the same as'
         print *, '         the number of orbitals listed in QRHF_GEN.'
         call errex
      end if

c ----------------------------------------------------------------------

c VALIDATE THE ASVs...

c   o print the ASV validation header
      if (bVerbose) then
         print '(/)'
         print *, '    ',
     &   '                  ACES STATE VARIABLE VALIDATION LOG'
         print *, '    ', ('-',iTmpReg11=1,70)
      end if

c   o set the default BASIS string
      if (ioppar(h_IOPPAR_basis).eq.0) BasNam='SPECIAL'

c   o identify the type of coordinates
      if (ioppar(h_IOPPAR_coordinates).eq.3) then
         if (bAutoCart) then
            szTmp = 'COORDINATES=CARTESIAN'//achar(0)
         else
            szTmp = 'COORDINATES=INTERNAL'//achar(0)
         end if
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c   o convert CALC aliases
      if (ioppar(h_IOPPAR_calc).eq.43.or.
     &    ioppar(h_IOPPAR_calc).eq.44    ) then
         if      (ioppar(h_IOPPAR_calc).eq.43) then
            call asv_update_kv('CALC=CCSD(T)'//achar(0),iASV,szData,0)
            call asv_update_kv('KUCHARSKI'//achar(0),iASV,szData,0)
         else if (ioppar(h_IOPPAR_calc).eq.44) then
            call asv_update_kv('CALC=CCSDT'//achar(0),iASV,szData,0)
            call asv_update_kv('KUCHARSKI'//achar(0),iASV,szData,0)
         else
            print *, '@GTFLGS: logic error at line ',__LINE__
            call errex
         end if
      end if

cJDW 9/16/97 Before too long, let's try to get a block of code which sets
c            reference and related keywords all in one place. Then we
c            can understand it better !
c
cMN MAKE RHF REFERENCE !!!
cMN      if (ioppar(h_IOPPAR_makerhf).eq.1) ioppar(h_IOPPAR_ref) = 0

c   o was QRHF_GENERAL set?
      if (ioppar(h_IOPPAR_qrhf_gen).ne.0) then
         call asv_update_kv('NON-HF'//achar(0),iASV,szData,0)
         DoQRHF = .true.
      else
         DoQRHF = .false.
      end if

      if (ioppar(h_IOPPAR_ref).ge.3) then
         call PutRec(20,'JOBARC','OSCALC  ',1,1)
         szTmp = 'REF=UHF'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         bTDCalc = .true.
      else
         bTDCalc = .false.
      end if

c NOW HAVE FINISHED READING IN THE STUFF.  MAKE SURE THAT MANDATORY
c  DEPENDENCIES ARE TAKEN CARE OF.  THIS BRANCH IS ALSO MADE IF THE
c  STRING IS READ IN THE OLD INTEGER LIST FORMAT.
c******************************************************************
c
c Ajith 03/2000, A bug fix to get the QRHF back on track!
c When ioppar(h_IOPPAR_makerhf) is turned on, Marcel would use RHF
c type structure (both Alpha and Beta have the same eigen
c vectors). However, if it is a normal QRHF calcualtion we
c need to force to be ROHF so that other NON-HF flags will
c get set properly. ioppar(h_IOPPAR_makerhf) is off by default. The
c other thing to worry about is the two determinant CC methods
c wich use QRHF key-word to specify the appropriate configurations.
c We can not run two determinant CC as ROHF since that would
c turn on the semi-canonical flag. The following logic will
c satisfy all the requirments and produce correct results for
c QRHF and two determinant calculations. Warning! The changes
c like these should be made with exterme caution and by only
c people who know what they are doing.

c If two determinant, turn it back to UHF. The orbital invariance
c of two determinant CC is a complicated issue. Please consult
c an expert if you have questions about these issues!
      if (DoQRHF) then
         if (ioppar(h_IOPPAR_makerhf).ne.0) then
            do i = 1, nIrMax
               nOcc(i,2) = nOcc(i,1)
            end do
         end if
         if (bTDCalc) then
            szTmp = 'REF=UHF'//achar(0)
         else
            szTmp = 'REF=ROHF'//achar(0)
         end if
         call asv_update_kv(szTmp,iASV,szData,0)
         iSpin = 2
      end if

      if (ioppar(h_IOPPAR_joda_print).gt.999) then
         szTmp = 'JODA_PRINT=999'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

      if (ioppar(h_IOPPAR_joda_print).lt.0) then
         szTmp = 'JODA_PRINT=0'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

      if (ioppar(h_IOPPAR_convergence).gt.10) then
         szTmp = 'CONV=10'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

      if (ioppar(h_IOPPAR_mult).lt.1) then
         print *, '@GTFLGS: Multiplicity must be >= 1'
         call errex
      end if
c
c   o when quadratic convergence is set we need to do 
c     RPP until the density difference fells below a 
c     certain threshold. This is how it is done without
c     introducing a new key-word (see vscf.F).
c     04/2006, Ajith Perera.
c
      if (ioppar(h_IOPPAR_scf_extrap) .eq. 2) then
         call iputrec(10, 'JOBARC', "RPP_1ST ", 1, 1)
      endif

c   o FINDIF_OLD? I thought this was killed. -Yau
      if (ioppar(h_IOPPAR_vib).eq.2.and.
     &    ioppar(h_IOPPAR_sym).gt.1     ) then
         call asv_update_kv('!sym'//achar(0),iASV,szData,0)
      end if

c NOW TAKE CARE OF A FEW SENSIBLE DEPENDENCIES. THESE ARE
c DONE INTERNALLY IF THE OPTCTL STRING IS READ IN.
c HERE, TURN CURVILINEAR TRANSFORMATION ON IF TRANSITION SEARCH IS BEING
c DONE, OR IF HESSIAN IS AVAILABLE AND VIBRATIONAL CALCULATION
c NOT REQUESTED.
      call GFName('FCM',FName,iLength)
      inquire(file=FName(1:iLength),exist=bTmp)
      bTmp = (bTmp                             .and.
     &        ioppar(h_IOPPAR_vib        ).eq.0.and.
     &        ioppar(h_IOPPAR_curvilinear).eq.0     )
      if (bTmp                             .or.
     &    ioppar(h_IOPPAR_opt_method ).eq.4.and.
     &    ioppar(h_IOPPAR_curvilinear).eq.0     ) then
         call asv_update_kv('curvilinear'//achar(0),iASV,szData,0)
c         print *, '@GTFLGS: Hessian will be transformed to ',
c     &            'curvilinear coordinates.'
      end if

c   o Kohn-Sham implies FOCK=AO, NON-HF=ON, ORBITALS=SEMICANONICAL
      if (ioppar(h_IOPPAR_ks_pot).ne.0) then
         szTmp = 'fock=ao'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         szTmp = 'orbitals=semicanonical'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         call asv_update_kv('non-hf'//achar(0),iASV,szData,0)
      end if

c   o vmol cannot do direct integrals
      if (ioppar(h_IOPPAR_integrals).eq.1.and.
     &    ioppar(h_IOPPAR_direct   ).eq.1     ) then
         print *, '@GTFLGS: vmol cannot calculate integrals directly'
         call errex
      end if

c   o direct integrals are only useful for FOCK=AO
      if (ioppar(h_IOPPAR_direct).eq.1) then
         szTmp = 'fock=ao'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c   o UNO ref implies semi-canonical orbitals, NON-HF=ON
      if (ioppar(h_IOPPAR_uno_ref).eq.1) then
         szTmp = 'orbitals=semicanonical'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         call asv_update_kv('non-hf'//achar(0),iASV,szData,0)
      end if

c   o put the UHF stuff in place
      if (ioppar(h_IOPPAR_ref).eq.0) then
         iUHF = 0
      else
         iUHF = 1
      end if
      call iputrec(20,'JOBARC','UHFRHF  ',1,iUHF)
      if (iUHF.eq.1.and.nDrop(1).ne.0) then
         nDrop(2) = nDrop(1)
         do i = 1, nDrop(1)
            iDrop(i,2) = iDrop(i,1)
         end do
      end if

c *****************************************
c ** INITIAL GEOMETRY OPTIMIZATION LOGIC **
c *****************************************

      if (bStar) then
c      o some coordinates have been flagged for relaxation
         if (ioppar(h_IOPPAR_geom_opt).eq.0) then
c         o elevate geom_opt from none to partial
            szTmp = 'geom_opt=partial'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      else
c      o no coordinates have been explicitly flagged for relaxation
         if (ioppar(h_IOPPAR_geom_opt).eq.0) then
c         o geom_opt=none
            if (ioppar(h_IOPPAR_opt_method).ne.0) then
               print *,
     &'@GTFLGS: A geometry relaxation algorithm has been selected,'
               print *,
     &'         but no asterisks appear in the Z-matrix and geom_opt is'
               print *,
     &'         not set to FULL.'
               call errex
            end if
         end if
         if (ioppar(h_IOPPAR_geom_opt).eq.1) then
c         o geom_opt=partial
            print *,
     &'@GTFLGS: A partial geometry optimization has been requested, but'
            print *,
     &'         but no asterisks appear in the Z-matrix.'
            if (ioppar(h_IOPPAR_coordinates).eq.1) then
               print *,
     &'         (Note: Partial geometry optimizations are meaningless'
               print *,
     &'          for systems defined in XYZ coordinates.)'
            end if
            call errex
         end if
      end if
c
c    o Geometry optimizations with exact Hessians need to be told how
c      the Hessian is calcualted. If not specified set it to 1 
c      (at each update Hessian is recalculated!).
c
      if (ioppar(h_IOPPAR_geom_opt)     .gt. 0 .and. 
     &    ioppar(h_IOPPAR_init_hessian) .eq. 3 .and. 
     &    ioppar(h_IOPPAR_eval_hess)    .eq. -1) then
          szTmp = 'eval_hess=0'//achar(0)
          call asv_update_kv(szTmp,iASV,szData,0)
      endif
c
c    o Geometry optimizations with exact Hessians need to be told
c      told to use it. If not, lets do it too. 07/2006, Ajith Perera.
c 
      if (ioppar(h_IOPPAR_geom_opt)     .gt. 0 .and.
     &    ioppar(h_IOPPAR_eval_hess)    .ge. 0 .and.
     &    ioppar(h_IOPPAR_init_hessian) .eq. 0) then
          szTmp = 'init_hessian=3'//achar(0)
          call asv_update_kv(szTmp,iASV,szData,0)
      endif
c
c    o Check to see whether we have a geometry and if not set
c      the jobarc record HAVEGEOM to 0 indicate the geometry optimizations
c      can proceed.
c
      Call igetrec(0, 'JOBARC', 'HAVEGEOM', i_length, i_havegeom)
c#ifdef _DEBUG_LVL0
      Print*, "A geometry for frequency {0,1} present:", i_length 
c#endif 
      if (i_length .lt. 0) Call iputrec(20, 'JOBARC', 'HAVEGEOM', 
     &                                 1, 0) 
c
c    o If the record exsist and it is non zero then either
c      the geometry optimization is successfuly completed (or
c      could even be a frequency calculation without prior 
c      optimization. 
c
      bGeomOpt = (ioppar(h_IOPPAR_geom_opt).ne.0)
c
      Call igetrec(-1, 'JOBARC', 'HAVEGEOM', 1, i_havegeom)
      if (i_havegeom .gt. 0) Then 
          szTmp = 'geom_opt=none'//achar(0)
          call asv_update_kv(szTmp,iASV,szData,0)
          bGeomOpt = (ioppar(h_IOPPAR_geom_opt).ne.0) 
      endif
c
c    o vibrational frequency calculations using a used defined
c      input, then also we have a geometry!
c
      if (.not. bGeomOpt .and. (ioppar(h_IOPPAR_vib).gt.0))  
     &    call iputrec(10, 'JOBARC', 'HAVEGEOM', 1, 1)
c      
      if (bGeomOpt.and.ioppar(h_IOPPAR_opt_method).eq.0) then
c      o set the default optimization algorithm to MANR (Perera 12/2001)
         szTmp = 'opt_method=manr'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c *********************************************
c ** END INITIAL GEOMETRY OPTIMIZATION LOGIC **
c *********************************************
C
C Stop all sorts of reorientations when SYM=NONE except 
C for vibrational frequancy calcualtions. Only the center
C of mass and the principal axis coordinate system separte
C Hamiltonian to vibrational, rotational and traslational
C motions. In arbitrary coordinate systems we need to 
C project those motions from the Hessian, and until that
C is done do not turn on the noreori for SYM=NONE
C runs. Ajith Perera, 12/08.
C
#ifdef _HES_PROJECT_ON
      If (ioppar(h_IOPPAR_sym).eq.0 .and.
     &    ioppar(h_IOPPAR_vib).ge.1 .and. 
     &    ioppar(h_IOPPAR_noreori) .eq. 2 ) then
          szTmp = 'NOREORI=ON'//achar(0)
          call asv_update_kv(szTmp,iASV,szData,0)
      Endif
#endif

c Turn off symmetry if subgroup=C1.
      if (ioppar(h_IOPPAR_subgroup).eq.1) then
         call asv_update_kv('!sym'//achar(0),iASV,szData,0)
      end if
      if ((ioppar(h_IOPPAR_geom_opt) .eq. 2) .or.
     &    (ioppar(h_IOPPAR_geom_opt) .eq. 3)) then
         szTmp = 'SYMMETRY=NONE'//achar(0)
         if (ioppar(h_IOPPAR_sym).ne.0) then
             call asv_update_kv(szTmp,iASV,szData,0)
         Endif
C
C For geo. opt with FULL/RIC and SYM=NONE turn on the NOREORI
C if it is not given as an input. Ajith Perera, 12/08
C 
          if (ioppar(h_IOPPAR_noreori) .eq. 2 .and. .not.
     &        ioppar(h_IOPPAR_vib).ge.1 ) then
             szTmp = 'NOREORI=ON'//achar(0)
             call asv_update_kv(szTmp,iASV,szData,0)
          endif
      endif
C
c Ajith 11/98.
c Raman intensity calculations:
c    SCF or MBPT(2) -> PROPS to second-order
c    CCSD           -> PROPS to EOM_NLO
      if (ioppar(h_IOPPAR_raman).eq.1) then
         if (ioppar(h_IOPPAR_calc).eq.0.or.
     &       ioppar(h_IOPPAR_calc).eq.1) then
            szTmp = 'props=second_order'//achar(0)
         else
         if (ioppar(h_IOPPAR_calc).eq.10) then
            szTmp = 'props=eom_nlo'//achar(0)
         end if
         end if
         call asv_update_kv(szTmp,iASV,szData,0)
      end if
c
c SCF only analytical second order properties (polarizabilities)
c with PROPS=EOM_NLO is possible (see setmet.f in cphf). So,
c now this PROPS=EOM_NLO and CALC has the proper correspondence
c (note: earlier version one has to use PROPS=SECOND_ORDER for
c SCF and PROPS=EOM_NLO for CCSD). 03/2006. Ajith Perera. 
c
cMN
c    LOGIC FOR DEALING WITH ioppar(h_IOPPAR_program)
c    (DETERMINES ACES2 OR MN_A3 CALCULATION)

c In case xmrcc is not part of the ACES suite, then one cannot run MN_A3
c calculations. At this point we need to set ACES2=.true.
      ACES2 = .true.
      ACES2 = .false.

c OPEN_SHELL, DEA_TDA, DIP_TDA, STEOM, and IP_EOM
      bOpenShell = .not.(ioppar(h_IOPPAR_ref    ).eq.0.or.
     &                   ioppar(h_IOPPAR_makerhf).eq.1.or.
     &                   ioppar(h_IOPPAR_uno_ref).eq.1)
      bTDA = (ioppar(h_IOPPAR_dip_calc).eq.1.or.
     &        ioppar(h_IOPPAR_dea_calc).eq.1)
      bSTEOM = (ioppar(h_IOPPAR_excite  ).eq.9.or.
     &          ioppar(h_IOPPAR_dip_calc).ge.1.or.
     &          ioppar(h_IOPPAR_dea_calc).ge.1)
      bMN_A3 = (bSTEOM                        .or.
     &          bTDA                          .or.
     &          ioppar(h_IOPPAR_ip_calc).ge.1 .or.
     &          ioppar(h_IOPPAR_calc   ).eq.40)

c RHF EA-EOM
      if (ioppar(h_IOPPAR_ea_calc).ge.1.and.
     &    .not.bOpenShell                   ) bMN_A3=.true.

c only SCF, MBPT(2), LCCSD, LCCD, CCD, CCSD, and ACCSD
      if (.not.(ioppar(h_IOPPAR_calc).eq. 0.or.
     &          ioppar(h_IOPPAR_calc).eq. 1.or.
     &          ioppar(h_IOPPAR_calc).eq. 5.or.
     &          ioppar(h_IOPPAR_calc).eq. 6.or.
     &          ioppar(h_IOPPAR_calc).eq. 8.or.
     &          ioppar(h_IOPPAR_calc).eq.10.or.
     &          ioppar(h_IOPPAR_calc).eq.40)
     &   ) ACES2 = .true.

c closed-shell parent state
cSSS      if (.not.(ioppar(h_IOPPAR_ref    ).eq.0.or.
cSSS     &          ioppar(h_IOPPAR_makerhf).eq.1.or.
cSSS     &          ioppar(h_IOPPAR_uno_ref).eq.1)
cSSS     &   ) ACES2 = .true.
c
      if (bOpenShell) ACES2 = .true.

c no properties
      if (ioppar(h_IOPPAR_props).ne.0) ACES2 = .true.

c excite=?
      if (ioppar(h_IOPPAR_excite).eq.2.or.
     &    ioppar(h_IOPPAR_excite).eq.6.or.
     &    ioppar(h_IOPPAR_excite).eq.8
     &   ) ACES2 = .true.

c no brueckner orbitals
cSSS      if (ioppar(h_IOPPAR_brueckner).ne.0.and.
cSSS     &    nDrop(1).gt.0                       ) bMN_A3=.true.

c    CHECK IF ACES CAN DO THE CALCULATION AT ALL

cMN analytical gradients are allowed in certain cases (not for STEOM yet)
c   Later it is decided by default if gradients are calculated analytically.

c      if (bSTEOM.and.
c     &    (.not.(ioppar(h_IOPPAR_opt_method).ge.5.or.
c     &           (ioppar(h_IOPPAR_opt_method).eq.2.and.
c     &            ioppar(h_IOPPAR_grad_calc ).ne.0)
c     &          )
c     &    )
c     &   ) ACES2 = .true.

      if (ACES2.and.bMN_A3) then
         print *, '@GTFLGS: Some combination of features is not'
         print *, '         implemented in MN_A3 or ACES2.'
         print *, ' STEOM -> NO TRIPLES, NO OPEN-SHELL REFERENCE '
         print *, ' NO DROPED CORE + BRUECKNER + OPEN SHELL.'
         call errex
      end if

      if (ioppar(h_IOPPAR_program).eq.2.and.bMN_A3) then
         print *, '@GTFLGS: ACES2 calculation requested, but ',
     &            'MN_A3 is required.'
         call errex
      end if

      if (ioppar(h_IOPPAR_program).eq.3.and.ACES2) then
         print *, '@GTFLGS: MN_A3 calculation requested, but ',
     &            'ACES2 is required.'
c         print *, ' CAN ONLY HAPPEN WHEN DEVELOPING CODE'
         call errex
      end if

c    SET PROGRAM FLAG h_IOPPAR_program
      if (ioppar(h_IOPPAR_program).eq.0) then
         if (bMN_A3.and.ACES2) then
C
cYAU - This was already killed a few conditionals up.
         else
            if (bMN_A3) then
               szTmp = 'program=mn_a3'//achar(0)
            else
               szTmp = 'program=aces2'//achar(0)
            end if
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
c     else if ([program != 0]) then
      else
         if (ioppar(h_IOPPAR_program).eq.2) then
            if (bMN_A3) then
               szTmp = 'program=mn_a3'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
         else
            if (ACES2) then
               print *, '@GTFLGS: I hope you are creating new options ',
     &                  'in MN_A3!'
            end if
         end if
c     end if ([program == 0])
      end if
c
c There are significant changes which are relevent to opt/freq 
c option and also to accomodate the change deriv_level options 
c (from {ANALYTICAL, NUMERICAL, AUTO} to  {NONE, ANALYTICAL, NUMERICAL})
c 04/2006, Ajith Perera
c
      ACES2 = ioppar(h_IOPPAR_program).eq.2
      bTmp = bGeomOpt.or.(ioppar(h_IOPPAR_vib).ge.2) 
      bMN_A3 = .not.ACES2

      if (bMN_A3) then
c
C Write a JOBARC record to tell Marcel that this is a new version.
C So he can write a wrapper to convert all the changes of flags
C to original ones that he enjoy using. This way he is immune to
C whatever the changes that I make to joda flags.
C Ajith Perera, 12/08.
C MN: This could be done always. Taking this out of the if statement.
c
         iversion = 270
         Call iputrec(20, "JOBARC", "ACES2_ID", 1, iversion)
C
C If a geometry optimization or vibrational frequency calculation
C (bTmp=.true.), turn on the analytical gradient flag
C
      if (ioppar(h_IOPPAR_grad_calc).eq.0 .and. bTmp) Then
         szTmp = 'grad_calc=analytical'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      endif
C
c      o set certain flags associated with analytical gradients
         eom = .not.bSTEOM.and.
     &         (ioppar(h_IOPPAR_excite ).ge.3.or.
     &          ioppar(h_IOPPAR_ea_calc).ge.4.or.
     &          ioppar(h_IOPPAR_ip_calc).ge.4)

c      o gradients are not available for qrhf/uno-ref/brueckner
c        calculations.

         if (DoQRHF .or. ioppar(h_IOPPAR_uno_ref).eq.1.or.
     &       ioppar(h_IOPPAR_brueckner).eq.1) NO_AGRAD=.True.
C
C If Gradient flag is set (to analytical), then reset it to based on
C the availability of analytical gradients.
C Ajith Perera, 12/2008.
C
         if (ioppar(h_IOPPAR_grad_calc).gt.0) then
c         o set grad_calc if none (default to analytical)
            if (No_AGRAD) then
               szTmp = 'grad_calc=numerical'//achar(0)
            else
               szTmp = 'grad_calc=analytical'//achar(0)
            end if
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
C
         if (bTmp.and.(ioppar(h_IOPPAR_grad_calc).eq.1)) then
            if (eom.or.ioppar(h_IOPPAR_print).gt.1) then
               szTmp = 'vtran=full'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
            if (eom) then
               szTmp = 'estate_prop=expectation'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
            call asv_update_kv('abcdfull'//achar(0),iASV,szData,0)
         end if
c     else if (ACES2) then
      else
         If (bTmp) then
            if (ioppar(h_IOPPAR_grad_calc).eq.0) then
                szTmp = 'grad_calc=analytical'//achar(0)
                call asv_update_kv(szTmp,iASV,szData,0)
            end if
         endif 
c     end if (bMN_A3)
      end if

c   o SET DERIVATIVE LEVEL and how gradient are evaluated 
c     frequency and geometry optimizations and other instantces
c     where grad_calc is set, but not the derivative level.
c
      if (ioppar(h_IOPPAR_deriv_lev).eq.-1) then
          if (ioppar(h_IOPPAR_grad_calc).ne.0) then
              szTmp = 'deriv=first'//achar(0)
              call asv_update_kv(szTmp,iASV,szData,0)
          end if
          if (ioppar(h_IOPPAR_vib) .eq.1) then
              szTmp = 'deriv=second'//achar(0)
              call asv_update_kv(szTmp,iASV,szData,0)
              if (ioppar(h_IOPPAR_grad_calc).eq.0) then
                  szTmp = 'grad_calc=analytical'//achar(0)
                  call asv_update_kv(szTmp,iASV,szData,0)
              endif 
          else if (ioppar(h_IOPPAR_vib) .gt.1) then
              szTmp = 'deriv=first'//achar(0)
              call asv_update_kv(szTmp,iASV,szData,0)
              if (ioppar(h_IOPPAR_grad_calc).eq.0) then
                  szTmp = 'grad_calc=analytical'//achar(0)
                  call asv_update_kv(szTmp,iASV,szData,0)
              endif
          end if
          if (ioppar(h_IOPPAR_geom_opt).ne.0) Then
              if (ioppar(h_IOPPAR_grad_calc).eq.0) then
                  szTmp = 'grad_calc=analytical'//achar(0)
                  call asv_update_kv(szTmp,iASV,szData,0)
              endif
              iTmpReg13 = -1
              iRecal = ioppar(h_IOPPAR_eval_hess)
CSSS              if (iRecal.ne.0) iTmpReg13 = mod(NCycle-1,iRecal)
CSSS              if (iRecal.eq.0) iTmpReg13 = 0
              if (ioppar(h_IOPPAR_vib).eq.1.or.iTmpReg13.eq.0) then
                  szTmp = 'deriv=second'//achar(0)
                  call asv_update_kv(szTmp,iASV,szData,0)
              end if
          end if
          if (ioppar(h_IOPPAR_resraman).eq.1) then
              szTmp = 'deriv=first'//achar(0)
              call asv_update_kv(szTmp,iASV,szData,0)
          end if
CSSS          if (ioppar(h_IOPPAR_excite).gt.0) Then
CSSS              szTmp = 'deriv=first'//achar(0)
CSSS              call asv_update_kv(szTmp,iASV,szData,0)
CSSS           endif
      end if
C The HF-DFT gradient calculations need solving Z-vector equations.
C So elevate the derivative level to second derivatives. The HF
C DFT grads will be required for geo. opt, vib. freq and single 
C point grads. P. Verma and A. Perera. 08/2008. 
C
          if (ioppar(h_IOPPAR_scf_type).eq.2) Then
C
              if (ioppar(h_IOPPAR_geom_opt)  .ne.0  .or. 
     &            ioppar(h_IOPPAR_vib)       .ne.0  .or.
     &            ioppar(h_IOPPAR_grad_calc) .eq.1) Then
C
                  szTmp = 'deriv=second'//achar(0) 
                  call asv_update_kv(szTmp,iASV,szData,0)
                  call asv_update_kv('save_ints'//achar(0),
     &                                iASV,szData,0)
              EndIF
C
          Endif 
C
#ifdef _DEBUG_LVL0
      Print*, "The value of grad calc", ioppar(h_IOPPAR_grad_calc)
      Print*, "The Deriv. level, havegeom?, and vib? :",
     &     ioppar(h_IOPPAR_deriv_lev),  i_havegeom, 
     &     ioppar(h_IOPPAR_geom_opt), ioppar(h_IOPPAR_vib) 
#endif
c
      If (ioppar(h_IOPPAR_vib)      .eq. 1 .and. 
     &    ioppar(h_IOPPAR_geom_opt) .ne. 0 .and. 
     &    i_havegeom .eq. 0)  Then
              szTmp = 'deriv=first'//achar(0)  
              call asv_update_kv(szTmp,iASV,szData,0)
      Endif
#ifdef _DEBUG_LVL0
      Print*, "The derivative level @exit", ioppar(h_IOPPAR_deriv_lev) 
#endif

c    
c    o If we have not yet change the default deriavative level from
c      -1 to 0, let's do it now. Certain things below this point
c      depends on having derivative level set to zero. 
c
      If (ioppar(h_IOPPAR_deriv_lev).eq.-1) then
         szTmp = 'deriv=zero'//achar(0) 
         call asv_update_kv(szTmp,iASV,szData,0) 
      Endif
c
c   o HF stability check can only be done for
c     single point energy calculations (not vib. freq.
c     or geom. opt.)
      if (ioppar(h_IOPPAR_hfstab  ).ne.0.and.
     &    ioppar(h_IOPPAR_geom_opt).eq.0) then
c      o save ints all the time for post scf and HF stab
c         if (ioppar(h_IOPPAR_hfstab).eq.2) then
            call asv_update_kv('save_ints'//achar(0),iASV,szData,0)
c         end if
      end if

c   o Brueckner calculation
      if (ioppar(h_IOPPAR_brueckner).ne.0) then
         call asv_update_kv('save_ints'//achar(0),iASV,szData,0)
         call asv_update_kv('non-hf'//achar(0),iASV,szData,0)
      end if

cKJW 1/18/98
c turn save_ints on for all fno methods
c props need to be first order for cc densities
      if (ioppar(h_IOPPAR_fno_keep).gt.0) then
         call asv_update_kv('save_ints'//achar(0),iASV,szData,0)
         if (ioppar(h_IOPPAR_calc).gt.1) then
            szTmp = 'props=first_order'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      end if

c   o TRAP METHODS WHICH DON'T WORK YET
      if (ioppar(h_IOPPAR_ref).lt.2) then
         if (iEnAva1(ioppar(h_IOPPAR_calc)).eq.0) then
            print *, '@GTFLGS: CALC=',ioppar(h_IOPPAR_calc),
     &               ' not implemented for RHF/UHF references.'
            call errex
         end if
      else
         if (ioppar(h_IOPPAR_ref).eq.2) then
            if (iEnAva2(ioppar(h_IOPPAR_calc)).eq.0) then
               print *, '@GTFLGS: CALC=',ioppar(h_IOPPAR_calc),
     &                  ' not implemented for ROHF references.'
               call errex
            end if
         end if
      end if

      if (bTDCalc) then
         if (iEnAva3(ioppar(h_IOPPAR_calc)).eq.0) then
            print *, '@GTFLGS: CALC=',ioppar(h_IOPPAR_calc),
     &               ' not implemented for TWODET references.'
            call errex
         end if
      end if

c   o If CC calc, then set MAXCYC and turn on RLE with order equal to 5
      if (ioppar(h_IOPPAR_calc).ge.5) then
         if (ioppar(h_IOPPAR_cc_maxcyc).eq.0) then
            szTmp = 'cc_maxcyc=50'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
         if (ioppar(h_IOPPAR_cc_exporder).eq.0) then
            szTmp = 'cc_exporder=5'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
         if (ioppar(h_IOPPAR_tamp_sum).eq.0) then
            szTmp = 'tamp_sum=5'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
         if (ioppar(h_IOPPAR_ntop_tamp).eq.0) then
            szTmp = 'ntop_tamp=15'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      else
         szTmp = 'cc_extrap=DIIS'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c   o PROPERTY CALCULATION FOR CORRELATED CALCULATION
#ifdef _DEBUG_LVL0
      Print*, "First order props and grad flags @enter:",
     &         ioppar(h_IOPPAR_props),ioppar(h_IOPPAR_grad_calc) 
#endif
      if (ioppar(h_IOPPAR_props).ge.1) then
         szTmp = 'deriv=first'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         if (ioppar(h_IOPPAR_grad_calc) .eq. 0) then
             szTmp = 'grad_calc=analytical'//achar(0) 
             call asv_update_kv(szTmp,iASV,szData,0) 
         endif 
      end if
#ifdef _DEBUG_LVL0
      Print*, "First order props, deriv level and grad_calc @exit:", 
     &         ioppar(h_IOPPAR_props), ioppar(h_IOPPAR_deriv_lev),
     &         ioppar(h_IOPPAR_grad_calc)
#endif
     
cJDW 6/16/95
c Try to keep ioppar(h_IOPPAR_deriv_lev)=1 for PROP=J_FC, J_SD, J_SO
cJDW/MN 10/23/95
c ioppar(h_IOPPAR_props).ne.13 added; JSC_ALL.
      if (ioppar(h_IOPPAR_props).ge. 2.and.
     &    ioppar(h_IOPPAR_props).ne.11.and.
     &    ioppar(h_IOPPAR_props).ne. 8.and.
     &    ioppar(h_IOPPAR_props).ne. 9.and.
     &    ioppar(h_IOPPAR_props).ne.10.and.
     &    ioppar(h_IOPPAR_props).ne.13     ) then
         szTmp = 'deriv=second'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         if (ioppar(h_IOPPAR_calc).gt. 3.and.
     &       ioppar(h_IOPPAR_calc).ne. 8.and.
     &       ioppar(h_IOPPAR_calc).ne.10.and.
     &       ioppar(h_IOPPAR_calc).ne.23     ) then
            print *, '@GTFLGS: Second-order properties not ',
     &               'available for CALC=',ioppar(h_IOPPAR_calc)
            call errex
         end if
      end if

c If this is a correlated ROHF calculation, then orbitals=semicanonical
c unless turned off explicitly. For second derivatives, however, they
c are turned off.
c SG 1/8/96
c Also set semicanonical for other non-HF triples calculations.
      if (ioppar(h_IOPPAR_orbitals).eq.-1) then
         if ((ioppar(h_IOPPAR_nonhf).eq.1    ).and.
     &       (ioppar(h_IOPPAR_calc).eq. 4.or.
     &        ioppar(h_IOPPAR_calc).eq.13.or.
     &        ioppar(h_IOPPAR_calc).eq.22    )
     &      ) then
               szTmp = 'orbitals=semicanonical'//achar(0)
         end if
         if (ioppar(h_IOPPAR_ref).eq.2) then
            if (ioppar(h_IOPPAR_calc     ).ge.1.and.
     &          ioppar(h_IOPPAR_deriv_lev).le.1     ) then
               szTmp = 'orbitals=semicanonical'//achar(0)
            end if
            if (ioppar(h_IOPPAR_calc     ).ge.1.and.
     &          ioppar(h_IOPPAR_deriv_lev).ne.0     ) then
               szTmp = 'orbitals=standard'//achar(0)
            end if
            if (ioppar(h_IOPPAR_calc).eq.0) then
               szTmp = 'orbitals=standard'//achar(0)
            end if
         else
            szTmp = 'orbitals=standard'//achar(0)
         end if
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c AO-LADDER DEFAULT DEPENDS ON MACHINE ARCHITECTURE!
      if (ioppar(h_IOPPAR_abcdtype).eq.2) then
c ???
      end if

c TRAP GRADIENT OR HESSIAN CALCULATIONS WHICH WON'T WORK FOR SOME REFS.
      if (ioppar(h_IOPPAR_deriv_lev).gt.0) then

         if (ioppar(h_IOPPAR_ref).lt.2) then
c         o check RHF/UHF
            if (iGrAva1(ioppar(h_IOPPAR_calc)).eq.0.and.
     &          ioppar(h_IOPPAR_deriv_lev)    .eq.1.and.
     &          ioppar(h_IOPPAR_grad_calc)    .eq.1     ) then
               print *, '@GTFLGS: RHF/UHF gradient calculations not ',
     &                  'possible for CALC=',ioppar(h_IOPPAR_calc)
               close(unit=LUZ,status='KEEP')
               call errex
            end if
            if (iHsAva1(ioppar(h_IOPPAR_calc)).eq.0.and.
     &          ioppar(h_IOPPAR_deriv_lev)    .eq.2     ) then
               print *, '@GTFLGS: RHF/UHF Hessian calculations not ',
     &                  'possible for CALC=',ioppar(h_IOPPAR_calc)
               close(unit=LUZ,status='KEEP')
               call errex
            end if
         else
         if (ioppar(h_IOPPAR_ref).eq.2) then
c         o check ROHF
            if (iGrAva2(ioppar(h_IOPPAR_calc)).eq.0.and.
     &          ioppar(h_IOPPAR_deriv_lev)    .eq.1.and.
     &          ioppar(h_IOPPAR_grad_calc)    .eq.1     ) then
               print *, '@GTFLGS: ROHF gradient calculations not ',
     &                  'possible for CALC=',ioppar(h_IOPPAR_calc)
               close(unit=LUZ,status='KEEP')
               call errex
            end if
            if (iHsAva2(ioppar(h_IOPPAR_calc)).eq.0.and.
     &          ioppar(h_IOPPAR_deriv_lev)    .eq.2     ) then
               print *, '@GTFLGS: ROHF Hessian calculations not ',
     &                  'possible for CALC=',ioppar(h_IOPPAR_calc)
               close(unit=LUZ,status='KEEP')
               call errex
            end if
c        end if ([ROHF])
         end if
c        end if ([RHF/UHF])
         end if
         if (bTDCalc) then
            if (iGrAva2(ioppar(h_IOPPAR_calc)).eq.0.and.
     &          ioppar(h_IOPPAR_deriv_lev)    .eq.1.and.
     &          ioppar(h_IOPPAR_grad_calc)    .eq.1     ) then
               print *, '@GTFLGS: TWODET gradient calculations not ',
     &                  'possible for CALC=',ioppar(h_IOPPAR_calc)
               close(unit=LUZ,status='KEEP')
               call errex
            end if
         end if

c        FOR GRADIENT CALCULATIONS INVOLVING TRIPLE EXCITATIONS,
c        USE ALWAYS PERTURBED CANONICAL OR SEMI-CANONICAL ORBITALS.
         if (ioppar(h_IOPPAR_pert_orb).eq.2) then
c         o reset UNKNOWN
            if ( ioppar(h_IOPPAR_calc).eq. 4.or.
     &           ioppar(h_IOPPAR_calc).eq. 9.or.
     &          (ioppar(h_IOPPAR_calc).ge.11.and.
     &           ioppar(h_IOPPAR_calc).le.22     )
     &         ) then
               szTmp = 'pert_orb=canonical'//achar(0)
            else
               szTmp = 'pert_orb=standard'//achar(0)
            end if
            call asv_update_kv(szTmp,iASV,szData,0)
         end if

c        Stop calculations for UHF/ROHF if triples and PERT_ORB=STANDARD.
c        This limitation will be removed eventually.
c        (By whom? -YAU)
         if (ioppar(h_IOPPAR_ref     ).ge.1.and.
     &       ioppar(h_IOPPAR_pert_orb).eq.0     ) then
            if ( ioppar(h_IOPPAR_calc).eq. 4.or.
     &           ioppar(h_IOPPAR_calc).eq. 9.or.
     &          (ioppar(h_IOPPAR_calc).ge.11.and.
     &           ioppar(h_IOPPAR_calc).le.22)
     &         ) then
               print *, '@GTFLGS: PERT_ORB=STANDARD for specific CALC'
               print *, '         or REFERENCE is not implemented.'
               call errex
            end if
         end if

c     end if ([deriv>0])
      end if

c   o RELAXED DENSITY -> NATURAL ORBITALS BY DEFAULT
c YAU - As far as I can grep, RDO is not used anywhere (and _relax_dens is dead)
c      if (ioppar(h_IOPPAR_rdo).eq.-1) then
c         if (ioppar(h_IOPPAR_relax_dens).eq.1) then
c            call asv_update_kv('rdo'//achar(0),iASV,szData,0)
c         else
c            call asv_update_kv('!rdo'//achar(0),iASV,szData,0)
c         end if
c      end if

c LOGIC FOR EOM-CC CALCULATIONS AND OTHER STUFF REQUIRING FORMATION OF H-BAR
c
cJDW 6/16/95
c Three lines for PROP=J_FC, J_SD, J_SO.
cMN/JDW 10/23/95
c Extra options for h_IOPPAR_props, h_IOPPAR_excite,
c h_IOPPAR_ea_calc, h_IOPPAR_ip_calc included.
      if (ioppar(h_IOPPAR_excite  ).eq.  3.or.
     &    ioppar(h_IOPPAR_excite  ).eq.  7.or.
     &    ioppar(h_IOPPAR_props   ).eq.  8.or.
     &    ioppar(h_IOPPAR_props   ).eq.  9.or.
     &    ioppar(h_IOPPAR_props   ).eq. 10.or.
     &    ioppar(h_IOPPAR_props   ).eq. 11.or.
     &    ioppar(h_IOPPAR_props   ).eq. 13.or.
     &    ioppar(h_IOPPAR_props   ).ge.100.or.
     &    ioppar(h_IOPPAR_ea_calc ).ge.  4.or.
     &    ioppar(h_IOPPAR_ip_calc ).ge.  4.or.
     &    ioppar(h_IOPPAR_dip_calc).ge.  2.or.
     &    ioppar(h_IOPPAR_dea_calc).ge.  2
     &   ) then

         if (ioppar(h_IOPPAR_estate_prop).eq.3) then
            print *, '@GTFLGS: Analytic response properties not ',
     &               ' available for EOM-CC calculations.'
            call errex
         end if

         if (ioppar(h_IOPPAR_calc).eq.0) then
            szTmp = 'deriv=second'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if

         if      (ioppar(h_IOPPAR_calc).eq. 0) then
            szTmp = 'eomref=none'//achar(0)
         else if (ioppar(h_IOPPAR_calc).eq. 1) then
            szTmp = 'eomref=mbpt(2)'//achar(0)
         else if (ioppar(h_IOPPAR_calc).eq.10.or.
     &            ioppar(h_IOPPAR_calc).eq.13.or.
     &            ioppar(h_IOPPAR_calc).eq.14.or.
     &            ioppar(h_IOPPAR_calc).eq.15.or.
     &            ioppar(h_IOPPAR_calc).eq.16.or.
     &            ioppar(h_IOPPAR_calc).eq.18.or.
     &            ioppar(h_IOPPAR_calc).eq.22.or.
     &            ioppar(h_IOPPAR_calc).eq.33.or.
     &            ioppar(h_IOPPAR_calc).eq.34.or.
     &            ioppar(h_IOPPAR_calc).eq.40
     &           ) then
            szTmp = 'eomref=ccsd'//achar(0)
         else
            print *, '@GTFLGS: EOM calculation not possible with ',
     &               'CALC=',ioppar(h_IOPPAR_calc)
            call errex
         end if
         call asv_update_kv(szTmp,iASV,szData,0)

         call asv_update_kv('hbar'//achar(0),iASV,szData,0)

cMN/JDW 8 EA
         if (ioppar(h_IOPPAR_ea_calc).ge.1) then
            call iputrec(20,'JOBARC','EASYM_A ',8,EARoot(1,1))
            if (ioppar(h_IOPPAR_ref).ge.1) then
               call iputrec(20,'JOBARC','EASYM_B ',8,EARoot(1,2))
            end if
         end if

cMN IP
         if (ioppar(h_IOPPAR_ip_calc).ge.1) then
            call iputrec(20,'JOBARC','IPSYM_A ',8,IPRoot(1,1))
            if (ioppar(h_IOPPAR_ref).ge.1) then
               call iputrec(20,'JOBARC','IPSYM_B ',8,IPRoot(1,2))
            end if
         end if
      end if

cMN
c New parameters in vee (does not work with old vee)
      if (ioppar(h_IOPPAR_zeta_type).eq.1) then
         if (ioppar(h_IOPPAR_zeta_conv).eq.0) then
            szTmp = 'zeta_conv=7'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
         if (ioppar(h_IOPPAR_zeta_maxcyc).eq.0) then
            szTmp = 'zeta_maxcyc=50'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      else
         if (ioppar(h_IOPPAR_zeta_conv).eq.0) then
            szTmp = 'zeta_conv=14'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      end if

      if (ioppar(h_IOPPAR_excite   ).ge.1.and.
     &    ioppar(h_IOPPAR_grad_calc).eq.1     ) then
         if (ioppar(h_IOPPAR_geom_opt ).ne.0.or.
     &       ioppar(h_IOPPAR_vib      ).ne.0.or.
     &       ioppar(h_IOPPAR_resraman ).eq.1.or.
     &       ioppar(h_IOPPAR_deriv_lev).eq.1
     &      ) then
            szTmp = 'estate_prop=unrelaxed'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      end if

c   o CALC=SCF FOR TDA AND =MBPT(2) FOR CIS(D)
      if (ioppar(h_IOPPAR_excite).eq.5.or.
     &    ioppar(h_IOPPAR_excite).eq.1    ) then
         szTmp = 'excite=tda'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         szTmp = 'calc=scf'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      else
         if (ioppar(h_IOPPAR_excite).eq.6) then
            szTmp = 'calc=mbpt(2)'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      end if

c VTRAN=FULL FOR CIS(D) GRADIENT CALCULATIONS
      if ((ioppar(h_IOPPAR_excite  ).eq.6.and.
     &     ioppar(h_IOPPAR_geom_opt).ne.0     ) .or.
     &    (ioppar(h_IOPPAR_excite  ).eq.6.and.
     &     ioppar(h_IOPPAR_vib     ).ne.0     )
     &   ) then
         szTmp = 'vtran=full'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

cMN set vtran/abcdfull to full for dea calculations.
      if (ioppar(h_IOPPAR_dea_calc).ne.0) then
         szTmp = 'vtran=full'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         call asv_update_kv('abcdfull'//achar(0),iASV,szData,0)
      end if

c SG 3/11/96
c  make sure that analytic gradients are available if requested
      if (ioppar(h_IOPPAR_excite   ).gt.0.and.
     &    ioppar(h_IOPPAR_deriv_lev).gt.0     ) then
         if (ioppar(h_IOPPAR_deriv_lev).gt.1) then
            print *, '@GTFLGS: Analytical frequencies are not',
     &               ' available for excited states.'
            call errex
         end if
         if (.not.( ioppar(h_IOPPAR_excite).eq.1.or.
     &              ioppar(h_IOPPAR_excite).eq.3.or.
     &              ioppar(h_IOPPAR_excite).eq.5.or.
     &              ioppar(h_IOPPAR_excite).eq.9.or.
     &             (ioppar(h_IOPPAR_excite).eq.7.and.
     &              ioppar(h_IOPPAR_eomref).eq.2)
     &            )
     &      ) then
            print *, '@GTFLGS: Analytic gradients are not available',
     &               ' for this excited state method.'
            call errex
         end if
      end if

c SG 5/7/96
c For excited state gradients and frequencies, set the convergence at
c   10**-7, otherwise set it at 10**-5
      if (ioppar(h_IOPPAR_estate_tol).eq.-1) then
         if ( ioppar(h_IOPPAR_excite  ).gt.0.and.
     &       (ioppar(h_IOPPAR_geom_opt).ne.0.or.
     &        ioppar(h_IOPPAR_vib     ).ne.0)
     &      ) then
            szTmp = 'estate_tol=7'//achar(0)
         else
            szTmp = 'estate_tol=5'//achar(0)
         end if
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

cYAU - unused
c 7002 FORMAT(T3,'@GTFLGS: Property calculations not compatible ',
c     &     'with analytic TDA gradients.')

c FOR HBAR CALCULATIONS, SET DERIVATIVE LEVEL TO 1
c      if (ioppar(h_IOPPAR_hbar     ).gt.0.and.
c          ioppar(h_IOPPAR_deriv_lev).eq.0     ) then
c         szTmp = 'deriv_lev=1'//achar(0)
c         call asv_update_kv(szTmp,iASV,szData,0)
c      end if

c SG 1/8/97 Set HF2_FILE=SAVE if ABCDTYPE=MULTIPASS
      if (ioppar(h_IOPPAR_abcdtype).eq.1) then
         szTmp = 'hf2_file=save'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

cMN DETERMINE INTEGRAL TRANSFORMATION
      if (ioppar(h_IOPPAR_eomref).ge.1.and.
     &    ioppar(h_IOPPAR_vtran ).eq.0     ) then
         szTmp = 'vtran=full'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

cYAU - This needs to be checked. The original IOPPAR index was called DIRECT,
c      but it was renamed to TURBOMOLE. It makes sense to switch to AOBASIS
c      for direct integrals, but I know nothing about how TURBOMOLE works so
c      I cannot say whether this is an actual dependency.
c SET ABCDTYPE=AOBASIS IF TURBOMOLE=ON
      if (ioppar(h_IOPPAR_turbomole).eq.1) then
         szTmp = 'abcdtype=aobasis'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c IF ABCDTYPE=AOBASIS, SET VTRAN=PARTIAL AND VICE VERSA
c ALSO SET GAMMA_ABCD=DIRECT
      if (ioppar(h_IOPPAR_vtran).eq.2) then
         szTmp = 'abcdtype=aobasis'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if
      if (ioppar(h_IOPPAR_abcdtype).eq.2) then
         szTmp = 'vtran=partial'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
         szTmp = 'gamma_abcd=direct'//achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c IF ROHF AND A CORRELATED CALCULATION, THEN SET NON-HF FLAG.
      if (ioppar(h_IOPPAR_ref).eq.2.and.ioppar(h_IOPPAR_calc).ge.1) then
         call asv_update_kv('non-hf'//achar(0),iASV,szData,0)
      end if

c HANDLE LOGIC FOR ENERGY-ONLY OPTIMIZATIONS AND FREQUENCY CALCULATIONS

c Ajith Perera 06/2001
c All the logic that deals with the FD_STEPSIZE is moved to
c the following block of code.
c      if (ioppar(h_IOPPAR_grad_calc).eq.1) then
      I2 = 0

c if (bGeomOpt) then; geometry optimization (no longer needed,
c here (see below), all we need to know whether we are doing a 
c finite difference gradient calculation. 
c
c   o set FD_STEPSIZE if the user has not
      if (bGeomOpt) then
c
         if (ioppar(h_IOPPAR_grad_calc).eq.2) then
c      o numerical gradients
            if (ioppar(h_IOPPAR_fd_stepsize).eq.0) then
               szTmp = 'fd_stepsize=25'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
         end if
c
         if (ioppar(h_IOPPAR_hess_update).eq.0) then
            if (ioppar(h_IOPPAR_opt_method).le.3) then
               ioppar(h_IOPPAR_hess_update) = 2
            else
               ioppar(h_IOPPAR_hess_update) = 4
            end if
         end if
c
c      o only symmetric distortions
         call iputrec(20,'JOBARC','NFDIRREP',1,1)
         call iputrec(20,'JOBARC','FDIRREP ',1,1)
c
      else if (ioppar(h_IOPPAR_vib).ge.2) then 
c      o vibrational frequency calculation
         if (ioppar(h_IOPPAR_grad_calc).eq.2) then
c          o numerical gradients
            if (ioppar(h_IOPPAR_fd_stepsize).eq.0) then
               szTmp = 'fd_stepsize=200'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
         else
            if (ioppar(h_IOPPAR_fd_stepsize).eq.0) then
               szTmp = 'fd_stepsize=50'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
         endif
c
      else if (ioppar(h_IOPPAR_grad_calc).eq.2) then
c      o numerical gradients 
         if (ioppar(h_IOPPAR_fd_stepsize).eq.0) then
            szTmp = 'fd_stepsize=25'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
c     end if (bGeomOpt)
      end if

c   o restarting something other than finite differences?
      if (ioppar(h_IOPPAR_restart).ne.0) then
         if (ioppar(h_IOPPAR_vib     ).ne.3.and.
     &       ioppar(h_IOPPAR_geom_opt).eq.0     ) then
            szTmp = 'restart=0'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      end if
c
cJDW 10/23/95. Block of code from MN.

c   LOGIC FOR DEALING WITH NEWVRT
      if (ioppar(h_IOPPAR_newvrt).eq.1.or.DoQRHF) then
cMN         if (ioppar(h_IOPPAR_orbitals).eq.1) then
cMN            print *, '@GTFLGS: SEMICANONICAL ORBITALS AND NEWVRT ',
cMN     &               'OPTIONS ARE INCOMPATIBLE'
cMN            call errex
cMN         end if
         call asv_update_kv('non-hf'//achar(0),iASV,szData,0)
      end if

cMN  LOGIC FOR DEALING WITH HBARABCD, HBARABCI
c
cJDW 10/23/95. Note that these settings are biased to the new vee code.
c              for the old vee code, we need DoABCD and DoABCI to be true
c              in the UHF/MO basis code.
c
c  SET HBARABCD/HBARABCI FLAGS: ONLY USED IN ACES2
c
c  FULL HBARABCD/HBARABCI IS REQUIRED FOR EA-EOMCC AND MOST VCCEH CALCULATIONS

      if (ACES2) then
c      o set DoABCD and DoABCI
         if ((ioppar(h_IOPPAR_props).ge. 8.and.
     &        ioppar(h_IOPPAR_props).le.11).or.
     &        ioppar(h_IOPPAR_props).eq.13
     &      ) then
            DoABCD = (ioppar(h_IOPPAR_excite ).lt.7.or.
     &                ioppar(h_IOPPAR_eomprop).gt.1)
            DoABCI = .true.
         else
            DoABCD = .false.
            DoABCI = .false.
         end if
c      o elevate DoABCD and DoABCI
         if (ioppar(h_IOPPAR_ea_calc).ge.5) then
            DoABCD = .true.
            DoABCI = .true.
         end if
         if (ioppar(h_IOPPAR_ea_calc).eq.4.and.
     &       ioppar(h_IOPPAR_eomref ).eq.1)
     &      DoABCI = .true.
         if (ioppar(h_IOPPAR_excite).eq.7.or.
     &       ioppar(h_IOPPAR_excite).eq.8)
     &      DoABCI = .true.
         if (ioppar(h_IOPPAR_eomref).eq.2.and.
     &       ioppar(h_IOPPAR_nonhf ).eq.0)
     &      DoABCI = .false.
c      o update ASVs
         if (DoABCD) then
            if (ioppar(h_IOPPAR_hbarabcd).eq.0) then
               call asv_update_kv('hbarabcd'//achar(0),iASV,szData,0)
            else
               if (ioppar(h_IOPPAR_hbarabcd).eq.1) then
cMN                  print *, '@GTFLGS: HBARABCD should be ON.'
               end if
            end if
         else
            if (ioppar(h_IOPPAR_hbarabcd).eq.0) then
               call asv_update_kv('!hbarabcd'//achar(0),iASV,szData,0)
            else
               if (ioppar(h_IOPPAR_hbarabcd).eq.2) then
cMN                  print *, '@GTFLGS: HBARABCD should be OFF.'
               end if
            end if
         end if
         if (DoABCI) then
            if (ioppar(h_IOPPAR_hbarabci).eq.0) then
               call asv_update_kv('hbarabci'//achar(0),iASV,szData,0)
            else
               if (ioppar(h_IOPPAR_hbarabci).eq.1) then
cMN                  print *, '@GTFLGS: HBARABCI should be ON.'
               end if
            end if
         else
            if (ioppar(h_IOPPAR_hbarabci).eq.0) then
               call asv_update_kv('!hbarabci'//achar(0),iASV,szData,0)
            else
               if (ioppar(h_IOPPAR_hbarabci).eq.2) then
cMN                  print *, '@GTFLGS: HBARABCI should be OFF.'
               end if
            end if
         end if
c     else if (.not.ACES2) then
      else
         call asv_update_kv('!hbarabcd'//achar(0),iASV,szData,0)
         call asv_update_kv('!hbarabci'//achar(0),iASV,szData,0)
c     end if (ACES2)
      end if

c   o Assymetric CCSD(T) or Lambda-based CCSD(T) designated as
c     ACCSD(T) requires to solve for CCSD/Lambda followed by
c     non-iterative triples correction with Lambda.
      if (ioppar(h_IOPPAR_calc).eq.42) then
         if (ioppar(h_IOPPAR_deriv_lev).eq.-1) then
            szTmp = 'deriv=zero'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         endif 
         if (ioppar(h_IOPPAR_grad_calc).eq.1) then
c         o This is an ACCSD(T) analytic derivative calculation,
c           then we need hbar to be formed after lambda.
            if (ioppar(h_IOPPAR_hbarabcd).eq.0) then
               call asv_update_kv('hbarabcd'//achar(0),iASV,szData,0)
            end if
            if (ioppar(h_IOPPAR_hbarabci).eq.0) then
               call asv_update_kv('hbarabci'//achar(0),iASV,szData,0)
            end if
         end if
      end if

cMN COMPRESS
c
c DETERMINE IF COMPRESSED ABCD INTEGRALS ARE TO BE USED.
c
c UNCOMPRESSED ABCD INTEGRALS ARE AT PRESENT REQUIRED FOR
c    1. CCSDT VARIANTS (SEE BELOW).
c    2. ACES2 EA_EOMCC CALCULATIONS.
c    3. ENERGY GRADIENT CALCULATIONS.
c    4. IN CASE OF EOMCC CALCULATIONS THE SETTING IS COMPLICATED.
c        WE FOLLOW THE LOGIC GIVEN BY HBARABCD/HBARABCI
c    5. CERTAIN EOM PROPERTY CALCULATIONS USING VCCEH.
c    6. DEA CALCULATIONS IN MN_A3
c    7. ANALYTICAL GRADIENTS IN MN_A3
c    8. ACCSD CALCULATIONS
c
c  THIS FLAG IS SET IF IT HAS THE DEFAULT VALUE ('UNKNOWN')
c  IN OTHER CASES WE LEAVE THE FLAG AS IS BUT PRINT OUT A WARNING
c
cJDW. 5/13/93.  bCCSDT3 determines whether the calculation is
c                  CCSDT-3, CCSDT-4, or CCSDT. This affects list 233 in
c                  RHF calculations. This variable occurs in the
c                  main program and in DS16AB.
cJDW. 6/28/93.  CCSDT3 extended to include noniterative fifth-order
c                  triples in "CCSD+T*(CCSD)" (CC5SD(T)).
cJDW. 10/14/93. CCSDT3 extended to include other noniterative
c                  fifth-order calculations.

c   o set bCompress
      bCCSDT3 = ( ioppar(h_IOPPAR_calc).eq.12      .or.
     &           (ioppar(h_IOPPAR_calc).ge.26.and.
     &            ioppar(h_IOPPAR_calc).le.31     ).or.
     &           (ioppar(h_IOPPAR_calc).ge.16.and.
     &            ioppar(h_IOPPAR_calc).le.18     ).or.
     &            ioppar(h_IOPPAR_calc).eq.33      .or.
     &            ioppar(h_IOPPAR_calc).eq.34
     &          )
      bCompress = (.not.bCCSDT3)
c   o elevate bCompress
      if (ioppar(h_IOPPAR_ea_calc).ge.5.and.ACES2)
     &   bCompress = .false.
      if (ioppar(h_IOPPAR_dea_calc).ge.1)
     &   bCompress = .false.
      if (ioppar(h_IOPPAR_calc).eq.40)
     &   bCompress = .false.
      if ((ioppar(h_IOPPAR_deriv_lev ).ne.  0.or.
     &     ioppar(h_IOPPAR_props     ).ne.  0     ).and.
     &    (ioppar(h_IOPPAR_excite    ).eq.  0.and.
     &     ioppar(h_IOPPAR_props     ).lt.100     ).and.
     &    (ioppar(h_IOPPAR_excite    ).eq.  0.and.
     &     ioppar(h_IOPPAR_props     ).ne. 11)          )
     &   bCompress = .false.
      if (ioppar(h_IOPPAR_excite     ).ge.3.and.
     &    ioppar(h_IOPPAR_estate_prop).eq.2)
     &   bCompress = .false.
      if (ioppar(h_IOPPAR_hbarabcd).eq.2.or.
     &    ioppar(h_IOPPAR_hbarabci).eq.2)
     &   bCompress = .false.
      if (bMN_A3)
     &   bCompress = .false.
c   o update ABCDFULL
      if (ioppar(h_IOPPAR_abcdfull).eq.0) then
         if (bCompress) then
            call asv_update_kv('!abcdfull'//achar(0),iASV,szData,0)
         else
            call asv_update_kv('abcdfull'//achar(0),iASV,szData,0)
         end if
      else
         if ((.not.bCompress.and.ioppar(h_IOPPAR_abcdfull).eq.2)) then
            print *
            print *, '@GTFLGS: COMPRESSED ABCD INTEGRALS NOT SUPPORTED'
            print *, '         SWITCH ABCDFULL TO ON'
            print *
            call errex
         end if
      end if
cMN END COMPRESS

c CHECK COMPATIBILITY COMPRESSED ABCD/AO ALGORITHM AND HBARABCD/HBARABCI
      if (ioppar(h_IOPPAR_hbarabcd).eq.2.or.
     &    ioppar(h_IOPPAR_hbarabci).eq.2    ) then
         if (ioppar(h_IOPPAR_abcdfull).eq.2) then
            print *, '@GTFLGS: INCOMPATIBLE COMPRESS,HBARABCD,HBARABCI'
            print *, '         JODA NEEDS TO BE FIXED!'
            call asv_update_kv('abcdfull'//achar(0),iASV,szData,0)
         end if
         if (ioppar(h_IOPPAR_abcdtype).eq.2) then
            print *, '@GTFLGS: INCOMPATIBLE AOBASIS,HBARABCD,HBARABCI'
            print *, '         JODA NEEDS TO BE FIXED!'
            szTmp = 'abcdtype=standard'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
         end if
      end if

cMN END

cKB/JDW 10/26/95. Put in block of code to set GUESS=CORE,PERT_ORB=
c                 CANONICAL for dropped core gradient/property/findif
c                 calculations. Also, set ORBITAL=SEMICANONICAL for
c                 ROHF. Write variable NDRGEO to JOBARC. This is 0
c                 if no dropped core and 1 if dropped gradient/prop,
c                 etc.
cJDW 8/26/97.     Removed restriction to core guess for dropped core.
c                 Changes made to scf code.
      NDRGEO = 0
      if (ioppar(h_IOPPAR_geom_opt ).ne.0 .or.
     &    ioppar(h_IOPPAR_grad_calc).eq.1 .or.
     &    ioppar(h_IOPPAR_vib      ).eq.2 .or.
     &    ioppar(h_IOPPAR_vib      ).eq.3 .or.
     &    ioppar(h_IOPPAR_props    ).eq.1
     &   ) then
         if (nDrop(1).gt.0) then
            szTmp = 'pert_orb=canonical'//achar(0)
            call asv_update_kv(szTmp,iASV,szData,0)
            if (ioppar(h_IOPPAR_ref).eq.2) then
               szTmp = 'orbitals=semicanonical'//achar(0)
               call asv_update_kv(szTmp,iASV,szData,0)
            end if
            NDRGEO = 1
         end if
      end if
      call iputrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

c   o reduce MEM to a multiple of 128 Bytes (16 floats)
      iTmpReg12 = ioppar(h_IOPPAR_mem)
      iTmpReg12 = iTmpReg12 / iIntFp
      iTmpReg12 = iTmpReg12 - iAnd(iTmpReg12,15)
      iTmpReg12 = iTmpReg12 * iIntFp
      if (ioppar(h_IOPPAR_mem).ne.iTmpReg12) then
         write(szTmp,*) 'mem=',iTmpReg12,achar(0)
         call asv_update_kv(szTmp,iASV,szData,0)
      end if

c SG 11/15/97
c Call FigIO to calculate CACHE_RECS and FILE_RECSIZ based on MEMORY_SIZE
      iTmpReg15 = ioppar(h_IOPPAR_cache_recs)
      iTmpReg16 = ioppar(h_IOPPAR_file_recsiz)
      call FigIO(ioppar(h_IOPPAR_mem),iTmpReg15,iTmpReg16)
      write(szTmp,*) 'cache_recs=',iTmpReg15,achar(0)
      call asv_update_kv(szTmp,iASV,szData,0)
      write(szTmp,*) 'file_recsiz=',iTmpReg16,achar(0)
      call asv_update_kv(szTmp,iASV,szData,0)

c   o print the ASV validation footer
      if (bVerbose) print *, '    ', ('-',iTmpReg14=1,70)

c ----------------------------------------------------------------------

c WARN THE USER OF ANY FINAL ASV COMBINATIONS OR JUST DIE

c If triples are requested, and standard orbitals are enforced by the user,
c then write out a warning message
      if (ioppar(h_IOPPAR_orbitals).eq.0.and.
     &    ioppar(h_IOPPAR_ref     ).eq.2.and.
     &    (ioppar(h_IOPPAR_calc).eq. 4.or.
     &     ioppar(h_IOPPAR_calc).eq.13.or.
     &     ioppar(h_IOPPAR_calc).eq.22    )
     &   ) then
         print *,
     &'\n',
     &'********** WARNING !!! *******************************\n',
     &'   The requested method is not correct to fourth order\n',
     &'nor satisfies the usual invariance of CC or MBPT\n',
     &'methods with respect to orbital rotations.\n',
     &'   We recommend using semicanonical orbitals!\n',
     &'********** WARNING !!! *******************************\n',
     &'\n'
      end if

cWJL 1/12/94 JDW; 3/14/94
c  Here we need to check that the user is doing a TDHF calculation only for
c  a RHF reference.  If not, then bomb out!
      if (ioppar(h_IOPPAR_ref ).ne.0.and.
     &    ioppar(h_IOPPAR_tdhf).eq.1     ) then
         print *, '@GTFLGS: TDHF calculation only valid for RHF.'
         call errex
      end if

cJDW 3/31/94
c Stop NMR calculations (apart from CCSDeH) which cannot use sphericals.
c Also, always stop if PROP=TDHF has been specified.
      if (ioppar(h_IOPPAR_spherical).eq. 1.and.
     &    ioppar(h_IOPPAR_props    ).ne. 7.and.
     &    ioppar(h_IOPPAR_props    ).ne.11.and.
     &    ioppar(h_IOPPAR_props    ).ge. 3     ) then
         print *, '@GTFLGS: SPHERICAL=ON impossible for this',
     &            ' kind of NMR calculation.'
         call errex
      end if
      if (ioppar(h_IOPPAR_props).eq.7) then
         print *, '@GTFLGS: Use TDHF=ON and $INPUT namelist for',
     &            ' TDHF calculations (see manual).'
         call errex
      end if

c   o SYM=FULL
      if (ioppar(h_IOPPAR_sym).eq.3) then
         if (ioppar(h_IOPPAR_integrals).eq.1) then
            print *, '@GTFLGS: VMOL does not understand high symmetry'
            call errex
         end if
         if (ioppar(h_IOPPAR_integrals).eq.4) then
            print *, '@GTFLGS: seward does not understand high symmetry'
            call errex
         end if
      end if

c   o SYM!=OFF
      if (ioppar(h_IOPPAR_integrals).eq.5.and.
     &    ioppar(h_IOPPAR_sym      ).gt.1     ) then
         print *, '@GTFLGS: GAMESS cannot use symmetry.'
         call errex
      end if

c   o VibFreq AND GeomOpt?
c      if (ioppar(h_IOPPAR_vib     ).ne.0.and.
c     &    ioppar(h_IOPPAR_geom_opt).ne.0     ) then
c         print *,
c     &'@GTFLGS: Geometry optimizations and vibrational frequencies'
c         print *,
c     &'         are currently mutually exclusive.'
c         call errex
c      end if

c   o Response density matrix makes no sense for analytical gradients.
      if (ioppar(h_IOPPAR_density).eq.1.and.
     &    ioppar(h_IOPPAR_grad_calc).eq.1.and.
     &    (ioppar(h_IOPPAR_geom_opt).ne.0.or.ioppar(h_IOPPAR_vib).ne.0)
     &   ) then
         print *, '@GTFLGS: Geometry optimizations or Frequency'
         print *, '         calculations are not valid with the'
         print *, '         response density matrix.'
         call errex
      endif

c   o INTEGRALS=CADPAC?
      if (ioppar(h_IOPPAR_integrals).eq.3) then
         print *, '@GTFLGS: Integral program not available.'
         call errex
      end if

cJDW 9/16/97.
c I am assuming ioppar(h_IOPPAR_ref) will not change further down!
c Trap situation when an open-shell occupation has been
c specified and the reference is RHF. Stop the job in this
c case --- in the past we used to set REFERENCE=UHF.
      if (ioppar(h_IOPPAR_occ).eq.2.and.ioppar(h_IOPPAR_ref).eq.0) then
         print *, '@GTFLGS: Open-shell OCCUPATION is incompatible'
         print *, '         with RHF reference. Please specify an'
         print *, '         appropriate reference.               '
         call errex
      end if

c   o (REF == RHF) && (MULT != 1)
      if (ioppar(h_IOPPAR_ref).eq.0.and.ioppar(h_IOPPAR_mult).ne.1) then
         print *, '@GTFLGS: MULTIPLICITY must be 1 for RHF.'
         call errex
      end if

c   o DROPMO limitations
      if (nDrop(1)+nDrop(2).ne.0) then
         if (ioppar(h_IOPPAR_abcdtype ).ne.0.and.
     &       ioppar(h_IOPPAR_deriv_lev).gt.0    ) then
            print *, '@GTFLGS: No AO-basis gradients with dropped MOs.'
            call errex
         end if
      end if
C
c   o NOREORI and finite difference vibrational frequency calculations.
c
      if (ioppar(h_IOPPAR_vib) .gt. 1) then
         if (ioppar(h_IOPPAR_noreori).eq.1) then
             Print*, " '@GTFLGS: NOREORI=ON is not compatible with"
             Print*, "finite difference frequency calculations"
             call errex
          endif
      endif
c ----------------------------------------------------------------------

c IF REQUESTED, PRINT A NICE TABULAR LISTING OF THE CONTROL PARAMETERS.
      if (iPrt.ne.0) then

c   o dump the handles and integers
      call asv_dump

c   o print the ASV strings header
      print *, '                        ACES STATE VARIABLES (STRINGS)'
      print *, '         ', ('-',iTmpReg17=1,60)

c   o BASIS
      print '(10x,a,a)', 'BASIS = ',BasNam(1:linblnk(BasNam))

c   o OCCUPATION
      if (ioppar(h_IOPPAR_occ).eq.0) then
         print '(10x,a)', 'OCCUPATION = [ESTIMATED BY SCF]'
      else
         print '(10x,a,8i4)', 'OCCUPATION = A: ',(nOcc(j,1),j=1,nIrrps)
         if (ioppar(h_IOPPAR_occ).eq.2) then
         print '(10x,a,8i4)', '             B: ',(nOcc(j,2),j=1,nIrrps)
         end if
      end if

c   o IP_SYM
      if (ioppar(h_IOPPAR_ip_calc).ne.0) then
         call igetrec(20,'JOBARC','IPSYM_A ',nIrrps,iArr)
         print '(10x,a,8i4)', 'IP_SYM = A: ',(iArr(j),j=1,nIrrps)
         if (ioppar(h_IOPPAR_ref).ne.0) then
         call igetrec(20,'JOBARC','IPSYM_B ',nIrrps,iArr)
         print '(10x,a,8i4)', '         B: ',(iArr(j),j=1,nIrrps)
         end if
      end if

c   o EA_SYM
      if (ioppar(h_IOPPAR_ea_calc).ne.0) then
         call igetrec(20,'JOBARC','EASYM_A ',nIrrps,iArr)
         print '(10x,a,8i4)', 'EA_SYM = A: ',(iArr(j),j=1,nIrrps)
         if (ioppar(h_IOPPAR_ref).ne.0) then
         call igetrec(20,'JOBARC','EASYM_B ',nIrrps,iArr)
         print '(10x,a,8i4)', '         B: ',(iArr(j),j=1,nIrrps)
         end if
      end if

c   o DROPMO
      if (nDrop(1).ne.0) then
         iTmpReg19 = min(12,nDrop(1))
         print '(10x,a,12i4)', 'DROPMO = ',(iDrop(j,1),j=1,iTmpReg19)
         if (nDrop(1).gt.12) then
            do while (iTmpReg19.lt.nDrop(1))
               iTmpReg20 = min(12,nDrop(1)-iTmpReg19)
               print '(19x,12i4)', (iDrop(iTmpReg19+j,1),j=1,iTmpReg20)
               iTmpReg19 = iTmpReg19 + iTmpReg20
            end do
         end if
      end if

c   o FD_IRREPS (max 8)
      if (ioppar(h_IOPPAR_fd_irreps).ne.0) then
         iTmpReg21 = ioppar(h_IOPPAR_fd_irreps)
         call igetrec(20,'JOBARC','FDIRREP ',iTmpReg21,iArr)
         print '(10x,a,8i4)', 'FD_IRREPS = ',(iArr(j),j=1,iTmpReg21)
      end if

c   o ESTATE_SYM (max 8)
      if (ioppar(h_IOPPAR_estate_sym).ne.0) then
         iTmpReg28 = ioppar(h_IOPPAR_estate_sym)
         call igetrec(20,'JOBARC','EESYMINF',iTmpReg28,iArr)
         print '(10x,a,8i4)', 'ESTATE_SYM = ',(iArr(j),j=1,iTmpReg28)
      end if

c   o QRHF ASVs
      if (ioppar(h_IOPPAR_qrhf_gen).ne.0) then
         iTmpReg22 = ioppar(h_IOPPAR_qrhf_gen)
         call igetrec(20,'JOBARC','QRHFIRR ',iTmpReg22,iArr)
         iTmpReg23 = min(12,iTmpReg22)
         print '(10x,a,12i4)', 'QRHF_GEN = ',(iArr(j),j=1,iTmpReg23)
         if (iTmpReg22.gt.12) then
            do while (iTmpReg23.lt.iTmpReg22)
               iTmpReg24 = min(12,iTmpReg22-iTmpReg23)
               print '(21x,12i4)', (iArr(iTmpReg23+j),j=1,iTmpReg24)
               iTmpReg23 = iTmpReg23 + iTmpReg24
            end do
         end if
      end if
      if (ioppar(h_IOPPAR_qrhf_orb).ne.0) then
         iTmpReg25 = ioppar(h_IOPPAR_qrhf_orb)
         call igetrec(20,'JOBARC','QRHFLOC ',iTmpReg25,iArr)
         iTmpReg26 = min(12,iTmpReg25)
         print '(10x,a,12i4)', 'QRHF_ORB = ',(iArr(j),j=1,iTmpReg26)
         if (iTmpReg25.gt.12) then
            do while (iTmpReg26.lt.iTmpReg25)
               iTmpReg27 = min(12,iTmpReg25-iTmpReg26)
               print '(21x,12i4)', (iArr(iTmpReg26+j),j=1,iTmpReg27)
               iTmpReg26 = iTmpReg26 + iTmpReg27
            end do
         end if
      end if

c   o print the ASV strings footer
      print *, '         ', ('-',iTmpReg18=1,60)

c     end if (iPrt.ne.0)
      end if

c ----------------------------------------------------------------------

c DUMP IOPPAR AND THE OCCUPATION VECTOR TO JOBARC.

      if (bExist_JOBARC) then

c      o IFLAGS
         call iputrec(20,'JOBARC','IFLAGS  ',dim_iflags,ioppar)
         call iputrec(20,'JOBARC','IFLAGS2 ',dim_iflags2,
     &                                             ioppar(dim_iflags+1))

c      o OCCUPATION
         call iputrec(20,'JOBARC','OCCUPYA ',nIrMax,nOcc(1,1))
         if (ioppar(h_IOPPAR_ref).ne.0) then
            call iputrec(20,'JOBARC','OCCUPYB ',nIrMax,nOcc(1,2))
         end if

c      o DROPMO
         call iputrec(20,'JOBARC','NUMDROPA',1,nDrop(1))
         if (nDrop(1).ne.0) then
            call iputrec(20,'JOBARC','MODROPA ',nDrop(1),iDrop(1,1))
         end if
         if (ioppar(h_IOPPAR_ref).ne.0) then
            call iputrec(20,'JOBARC','NUMDROPB',1,nDrop(2))
            if (nDrop(2).ne.0) then
               call iputrec(20,'JOBARC','MODROPB ',nDrop(2),iDrop(1,2))
            end if
         end if

c      o JODAFLAG (17 consecutive ASVs)
         call iputrec(20,'JOBARC','JODAFLAG',JPARAM,
     &               ioppar(h_IOPPAR_joda_print))
      end if

      i1  = ioppar(h_IOPPAR_joda_print)
      i2  = ioppar(h_IOPPAR_opt_method)
      i3  = ioppar(h_IOPPAR_convergence)
      i4  = ioppar(h_IOPPAR_eigenvector)
      i5  = ioppar(h_IOPPAR_negeval)
      i6  = ioppar(h_IOPPAR_curvilinear)
      i7  = ioppar(h_IOPPAR_stp_siz_ctl)
      i8  = ioppar(h_IOPPAR_max_step)
      i9  = ioppar(h_IOPPAR_vibration)
      i10 = ioppar(h_IOPPAR_eval_hess)
      i11 = ioppar(h_IOPPAR_integrals)
      i12 = ioppar(h_IOPPAR_fd_stepsize)
      i13 = ioppar(h_IOPPAR_points)
      i14 = ioppar(h_IOPPAR_contraction)
      i15 = ioppar(h_IOPPAR_symmetry)
      i16 = ioppar(h_IOPPAR_basis)
      i17 = ioppar(h_IOPPAR_spherical)
cYAU - already initialized to zero
c      iErr = 0

c stop if RESET_FLAGS is ON
      if (ioppar(h_IOPPAR_reset_flags).eq.1) then
         call aces_ja_fin
         call joda_exit(0, '@GTFLGS')
      end if

c      close(unit=LUZ,status='KEEP')
      return

c ----------------------------------------------------------------------

 8000 write(*,*) '@GTFLGS: There was a problem reading ZMAT.'
cYAU - Why set this and then crash?
c      iErr = 1
      call errex

cYAU - not used
c 8111 iErr = 2
c      GOTO 159

 5400 write(*,*) '@GTFLGS: ZMAT is missing the ACES2 namelist.'
      close(unit=luz,status='KEEP')
      call errex

      return
      end

