C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine Driver(LUInt, NVMol, NAO, NSO, NCen, NOrbits, IOrdGp,
     $   SOAO, CenLbl, AngLbl, NAOOrb, NSOOrb, IAngBF, NAOWrk,
     $   IAOCen, NCenW1, NCenW2, Coord1, Coord2, CenMap, LScr, Scr,
     $   Scr2)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: driver.f,v 1.2 2005/03/09 16:50:21 yau Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     driver -- driver routine for VMol to JOBARC translation
C
C SYNOPSIS
      Integer LUInt, NAO, NSO, NCen, NOrbits, IOrdGp, Map
      Logical NVMol
      Integer LScr
      Double precision Scr(LScr, LScr), Scr2(LScr, LScr)
      Double precision SOAO(NAO, NAO)
      double precision xone,xzero
      Character*4 CenLbl(NSO), AngLbl(NSO)
      Integer NAOOrb(NSO), NSOOrb(NCen), IAngBF(NAO), NAOWrk(NAO)
      Integer IAOCen(NAO)
      Integer NCenW1(NCen), NCenW2(NCen), CenMap(NCen)
      Double precision Coord1(3, NCen), Coord2(3, IOrdGp-1)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Include 'mach.inc'
C $Id: driver.f,v 1.2 2005/03/09 16:50:21 yau Exp $
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C     Include 'flags.inc'
C $Id: driver.f,v 1.2 2005/03/09 16:50:21 yau Exp $
      Integer IFlags(100), JFlags(16)
      Common /Flags/ IFlags
      Common /Flags3/ JFlags
C LOCAL VARIABLES
C     I       General use
C     L       General use
C     CmpPgp  Computational point group
C     IOrdGp  Order of the computational point group
C     S       Scaling factor for inverting SOAO
C     SEWARD  If it is true SEWARD integrals are being used.
C     ACESIII  If it is true ACESIII is used.
C
      Integer I, L, IACESIII
      LOGICAL SEWARD, ACESIII
      Character*4 CmpPgp
      Double precision S
C
        iErr = 0
      SEWARD = IFLAGS(56).EQ.4
      CALL GETREC(-1, "JOBARC", "AIIIJARC", 1, IACESIII)
      IF (IACESIII .GT. 0) ACESIII = .TRUE.
C
C
C     ***********************************
C     * Get the SO -> AO transformation *
C     ***********************************
C     "New VMol" supports spherical harmonic basis sets.  These generally
c     have rectangular SOAO matrices (because "contaminants" are usually
C     deleted) which cause problems when the time comes to invert them.
c     Fortunately, New VMol happens to write the SOAO entries for the
C     deleted orbitals after the true SOAO matrix (because I made it :-).
C     We need the whole matrix to do this stuff, so if its a New VMol
C     run, we use the dimensions to insure the whole matrix is read.
C
C     "Old VMol" never writes anything other than a NAOxNAO matrix.
C
C     Write both the "full" SOAO matrix out, in case anyone needs it,
C     and the regular rectangular one.  The record AO2SO is a misnomer
C     and should be changed (to SO2AO)!
C
      IF(SEWARD .OR. ACESIII)THEN
       CALL GETREC(20,'JOBARC','FULLSOAO',NAO*NAO*IINTFP,SOAO)
       Call PutRec( 20, 'JOBARC', 'AO2SO   ', NAO*NSO*IIntFP, SOAO)
      ELSE
       If (NVMol) then
          Call RdSOAO(LUInt, 'SYMTRANS', NVMol, NAO, NAO, SOAO, NAO)
          Call PutRec( 20, 'JOBARC', 'AO2SO   ', NAO*NSO*IIntFP, SOAO)
       Else
          Call RdSOAO(LUInt, 'SYMTRANS', NVMol, NAO, NAO, SOAO, NAO)
          Call PutRec( 20, 'JOBARC', 'AO2SO   ', NAO*NAO*IIntFP, SOAO)
       EndIf
       Call PutRec( 20, 'JOBARC', 'FULLSOAO', NAO*NAO*IIntFP, SOAO)
      ENDIF
C
cSSS      write(6,*) ' @DRIVER-I, Matrix SOAO (SYMTRANS) '
cSSS      CALL OUTPUT(SOAO,1,NAO,1,NAO,NAO,NAO,1)
C
C     **********************************
C     * Invert the SOAO transformation *
C     **********************************
C     The full inverse can, of course act to the left or right of SOAO.
C     The submatrix can only act on the LEFT of SOAO!
C
C     The record name of the submatrix, AO2SOINV, has been chosen to
C     complement the (misnamed) AO2SO.  This should be changed (to 
C     AO2SO) when AO2SO is changed!
C
C     Copy the full matrix this time, and invert it.
C
      Call SCopy(NAO*NAO, SOAO, 1, Scr, 1)
c OLD
c      Call MInv(Scr, NAO, NAO, Scr2, S, 1.0d-30, 0, 1)
c NEW
      iErr = 0
      Call DGETRF(NAO,NAO,Scr,NAO,Scr2,iErr)
      if (iErr.ne.0) then
         if (iErr.gt.0) then
            print *, '@DRIVER: factorized a singular matrix'
            if (iflags(62).eq.1) then
            print *, '         Restart with SPHERICAL=OFF.'
            else
            print *, '         Restart with SYMMETRY=OFF.'
            end if
         else
            print *, '@DRIVER: dgetrf argument ',-iErr,' is illegal'
         end if
         call errex
      end if
      Call DGETRI(NAO,Scr,NAO,Scr2,Scr2(1,2),LScr*(LScr-1),iErr)
      if (iErr.ne.0) then
         if (iErr.gt.0) then
            print *, '@DRIVER: inverted a singular matrix'
            if (iflags(62).eq.1) then
            print *, '         Restart with SPHERICAL=OFF.'
            else
            print *, '         Restart with SYMMETRY=OFF.'
            end if
         else
            print *, '@DRIVER: dgetri argument ',-iErr,' is illegal'
         end if
         call errex
      end if
c END
C
C     Write the full inverse and the submatrix without the deleted
C     functions
C
      Call PutRec( 20, 'JOBARC', 'FULLAOSO', NAO*NAO*IIntFP, Scr)
      Call BlkCpy2(Scr, NAO, NAO, Scr2, NSO, NAO, 1, 1)
      Call PutRec( 20, 'JOBARC', 'AO2SOINV', NSO*NAO*IIntFP, Scr2)
C
C     ******************************************************
C     * Create computational <-> cartesian transformations *
C     ******************************************************
C     Only need to do this if the computational basis is _not_
C     the full cartesian basis -- this is the only time the
C     CSYMTRAN record will be written by VMOL.  Otherwise the
C     transformation is the unit matrix, and we don't bother to
C     write it out.
C
      If (IFlags(62) .eq. 1 .OR. NSO .ne. NAO) then
C
C        Get the full cartesian SO -> AO transformation
C
         IF(SEWARD .OR. ACESIII)THEN
          CALL GETREC(20,'JOBARC','CSYMTRAN',NAO*NAO*IINTFP,SOAO)
         ELSE
          Call RdSOAO(LuInt, 'CSYMTRAN', NVMol, NAO, NAO, SOAO, NAO)
         ENDIF
C
cSSS         write(6,*) ' @DRIVER-I, Matrix SOAO (CSYMTRAN) '
cSSS         CALL OUTPUT(SOAO,1,NAO,1,NAO,NAO,NAO,1)
C
C        AO -> computational SO transformation is in Scr from above.
C        Multiply by cartesian SO -> AO on the right to get the
C        cartesian SO (columns) -> computational SO (rows) 
C        transformation.
C
         xone=1.0d0
         xzero=0.0d0
         Call XGEMM('N', 'N', NAO, NAO, NAO, xone, Scr, NAO,
     $      SOAO, NAO, xzero, Scr2, NAO)
C
C        Write the transformation including orbitals dropped in the
C        computational basis (CART3CMP) as well as the submatrix
C        corresponding to the true computational basis (CART2CMP).
C
         Call PutRec( 20, 'JOBARC', 'CART3CMP', NAO*NAO*IIntFP, Scr2)
         Call BlkCpy2(Scr2, NAO, NAO, Scr, NSO, NAO, 1, 1)
         Call PutRec( 20, 'JOBARC', 'CART2CMP', NSO*NAO*IIntFP, Scr)
C
cSSS         Call DGEWRL('N', NAO, NAO, Scr2, NAO, 0, 0, 'junk', 'junk',
cSSS     $      'Num', 'Num', 1, 6, 'Computational to cartesian', 80,
cSSS     $      'F12.6', I)
C
C        Once again, a non-unitary transformation, so we must invert it.
C
C        Then write the full (CMP3CART) and true (CMP2CART)
c        computational basis results.
C
c OLD
c         Call MInv( Scr2, NAO, NAO, Scr, S, 1.0d-30, 0, 1)
c NEW
         Call DGETRF(NAO,NAO,Scr2,NAO,Scr,iErr)
         if (iErr.ne.0) then
            if (iErr.gt.0) then
               print *, '@DRIVER: factorized a singular matrix'
               print *, '         Restart with SPHERICAL=OFF.'
            else
               print *, '@DRIVER: dgetrf argument ',-iErr,' is illegal'
            end if
            call errex
         end if
         Call DGETRI(NAO,Scr2,NAO,Scr,Scr(1,2),LScr*(LScr-1),iErr)
         if (iErr.ne.0) then
            if (iErr.gt.0) then
               print *, '@DRIVER: inverted a singular matrix'
               print *, '         Restart with SPHERICAL=OFF.'
            else
               print *, '@DRIVER: dgetri argument ',-iErr,' is illegal'
            end if
            call errex
         end if
c END
         Call PutRec( 20, 'JOBARC', 'CMP3CART', NAO*NAO*IIntFP, Scr2)
         Call PutRec( 20, 'JOBARC', 'CMP2CART', NAO*NSO*IIntFP, Scr2)
C
cSSS         Call DGEWRL('N', NAO, NAO, Scr2, NAO, 0, 0, 'junk', 'junk',
cSSS     $      'Num', 'Num', 1, 6, 'Inverse computational to cartesian',
cSSS     $      80, 'F12.6', I)
C
C        Restore the computational SO -> AO matrix, needed below
C
         Call GetRec( 20, 'JOBARC', 'FULLSOAO', NAO*NAO*IIntFP, SOAO)
       EndIf
C 
C     **********************************
C     * Read the basis function labels* 
C     ********************************** 
C
      Call RdBLab(LUInt, NVMol, NSO, CenLbl, AngLbl)
C
C     *******************************
C     * Analyze the SO basis labels *
C     *******************************
C     Uses SOAO and the SO basis functions labels to evaluate
C     ANGMOMBF, NBASATOM, and NAOBFORB
C
      Call GetAOSO(NVMol, NAO, NSO, SOAO, CenLbl, AngLbl, NAOOrb, 
     $   IAngBF, NAOWrk, NCen, NOrbits, NSOOrb)
C
      Call PutRec( 20, 'JOBARC', 'ANGMOMBF', NAO, IAngBF)
      Call PutRec( 20, 'JOBARC', 'NBASATOM', NOrbits, NSOOrb)
      Call PutRec( 20, 'JOBARC', 'NAOBFORB', NOrbits, NAOOrb)
C
C     ****************************************************
C     * Determine the mapping from VMol to ZMAT ordering *
C     ****************************************************
C     Uses a bunch of info about the centers to figure out the
C     relation between VMOl's center ordering and the ZMAT input.
C
      L = 1
      Call GetCRec( -1, 'JOBARC', 'COMPPTGP', 4, CMPPGP)
      Call GetRec( -1, 'JOBARC', 'COMPPOPV', NOrbits, NCenW1)
      Call GetRec( -1, 'JOBARC', 'COMPMEMB', NCen, NCenW2)
      L = 3 * NCen * IIntFP
      Call GetRec( -1, 'JOBARC', 'COORD   ', L, Coord1)
C
      Call Remap( CMPPGP, IORDGP, NCen, NOrbits, Coord1, Coord2,
     &               NCenW2, NCenW1, CenMap)
C
      Call PutRec( 20, 'JOBARC', 'MAP2ZMAT', NCen, CenMap)
C
C
C     ************************************************
C     * Imitate VMol for mapping from AOs to centers *
C     ************************************************
C     Requires the population of each orbit in the computational basis
C
      Call VBFCen(NOrbits, NCenW1, NAOOrb, IAOCen)
      Call PutRec(20, 'JOBARC', 'CENTERBF', NAO, IAOCen)
cSSSC
cSSSC     Reorder CENTERBF, already in core.
cSSSC
cSSS      Call ISctr(NAO, IDummy, ICheck, IAngCmp)
cSSS      Call PutRec(20, 'JOBARC', 'CNTERBF0', NAO, IAngCmp)
cSSS      Write (6, 9000) 'CNTERBF0', (i, IAngCmp(i), i = 1, NAO)
cSSSC
cSSSC     Reorder ANGMOMBF
cSSSC
cSSS      Call GetRec( 20, 'JOBARC', 'ANGMOMBF', NAO, IDummy)
cSSS      Call ISctr(NAO, IDummy, ICheck, IAngCmp)
cSSS      Call PutRec(20, 'JOBARC', 'ANMOMBF0', NAO, IAngCmp)
cSSS      Write (6, 9000) 'ANMOMBF0', (i, IAngCmp(i), i = 1, NAO)
cSSSC
cSSSC     Reorder the transformation; written back out by the caller.
cSSSC
cSSS      Do 100 I = 1, NAO
cSSS         Call SCOPY(NBas, Evec(i, 1), NAO, Scr( ICheck(i), 1), NAO)
cSSS 100  Continue
cSSS      Call SCOPY(NAO*NBas, Scr, 1, Evec, 1)
cSSSC
cSSS 9000 Format(1X,'new style ', A,':'/(2I5))
      Return
      End
