C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine v2ja(Icore, Icrsiz)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: v2ja.F,v 1.6 2008/11/19 15:48:33 perera Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     v2ja -- write some VMol mapping information to CRAPS JOBARC file
C
C DESCRIPTION
C     A number of operations in CRAPS are dependent upon knowing the
C     relationship between the atomic centers and the basis set as
C     used by the integral program and as used by JODA and cohorts.
C     The program reads a VMol integral file to obtain the necessary
C     information and writes it on the JOBARC file for use by CRAPS
C     programs.  In this way, some of the more esoteric dependencies
C     of CRAPS of the chosen integral package can be removed from other
C     programs, where they don't really belong.
C
C     JOBARC records written by this program include:
C
C     NUCREP    VMOl nuclear repulsion energy
C     NUMBASIR  Nr. basis functions per irrep
C     NBASTOT   Total nr. of basis functions
C               The two records immediately above are numbers of
C               symmetry adapted orbitals -- the actual computational
C               basis.
C     NAOBASFN  Number of AO basis functions.  May be different from
C               NBASTOT if (for example) spherical harmonics are
C               used, since the raw AO basis is always in cartesians.
C
C     CENTERBF  Center number corresponding to each fn., VMol AO order
C     ANGMOMBF  l quantum number of each fn., VMOl AO order.
C     NBASATOM  Nr. SO basis fns. in each orbit.
C     NAOBFORB  Nr. AO basis fns in each orbit.
C
C     AO2SO     Misnamed record containing the SO -> AO transformation
c               matrix from VMol.  Not normalized.  May be rectangular.
C               To get the inverse of this matrix, divide each column
C               by the sum of the squares of the elements of the column.
C               The transpose of this matrix is the (left) inverse of AO2SO
C
C     MAP2ZMAT  Mapping function for atomic centers, VMol -> ZMAT (JODA)
C     CNTERBF0  CENTERBF in ZMAT ordering
C     ANMOMBF0  ANGMOMBF in ZMAT ordering
C
C     CMP2ZMAT  Transformation from VMol SO to ZMAT AO swapped ordering.
C               Multiplies on the left of the quantity being transformed.
C               NOT unitary.
C     ZMAT2CMP  Inverse of CMP2ZMAT.  NOT unitary.
C
C FUTURE IDEAS
C     normalized AO2SO matrix
C     direct product table
C     Labels for irreps
C     Cartesian->spherical transformation
C     m angular momentum labels
C
C ROUTINES REQUIRED
      External ErrEx, AOSize, CrapsI, PutRec, GetRec, IZero, aces_fin,
     $   InsMem, Driver, MkC2Z1, PDSwap, Transp, SNrm2, SScal, SDot
      Double precision SNrm2, SDot
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  NOTES ON SEWARD/ALASKA/MCKINLEY COMPATIBILITY. Ajith Perera, 07/2000
C
C vmol2ja has been modified to work with Seward/Alaska (Molcas). During
C Molcas runs, the integral program (SEWARD) does not create any
C SO->AO or AO->SO transformation matrices. Therefore, they are created
C here in GenTranMat. GenTranMat is a hacked up version of the old vmol
C Readin routine. The use of these transformation matrices in the rest
C of the code is sparse. They are listed here:
C  
C   AO2SO (misnomer) : transforms SYMMETRY-ADAPTED orbitals (Cartesian
C                      or Spherical) to NON-SYMMETRY-ADAPTED orbitals
C                      (Cartesian basis)
C                      Used in vscf, props, vee.
C   AO2SOINV         : the inverse transformation of AO2SO
C                      By default, it acts FROM THE LEFT.
C                      Used in props and vee. 
C
C   FULLSOAO         : transforms symmetry-adapted CARTESIAN AOs to
C                      non-symmetry-adapted CARTESIAN AOs
C                      Not used anywhere (Good). In Spherical computational
C                      basis set, the droped Cartesian functions are at the
C                      bottom of these matrices.
C   FULLAOSO         : the inverse transformation of FULLSOAO
C                      Also not used.
C
C The use of AO2SO and its inverse follows,
C
C   AO2SO(NAO,NSO)*MAT(NSOxNSO)*Trans(AO2SO) -> MAT(NAOxNAO))
C
C   AO2SOINV(NSO,NAO)*MAT(NAOxNAO)*Trans(AO2SOINV) -> MAT(NSOxNSO))
C
C             (NOTE THAT ALL THE QUANTITIES ARE IN VMOL ORDER)  
C             (NOTE THAT NAO REFERS TO THE NUMBER OF CARTESIAN BASIS FUNCTIONS)
C
C These four records are identical for Cartesian functions 
C and become the unit matrix when SYMMETRY=OFF.
C 
C   CART3CMP : transforms Cartesian symmetry adapted basis functions to
C              Cartesian symmetry adapted basis functions (correct)
C              For spherical computational basis sets, the deleted
C              Cartesian spherical harmonic contaminants (symmetry 
C              adapted) are stored from the last row up.
C              For Cartesian only cases, this is the unit matrix.
C              (not used, GOOD!)
C
C   CART2CMP : transforms a Cartesian symmetry adapted basis functions to
C              computational symmetry adapted basis functions. The
C              computational basis can be Cartesian or spherical depending on
C              users choice. (for Cartesian only cases this is a unit matrix)  
C              (not used)
C
C The inverse transformations are CMP3CART and CMP2CART. They do exactly the
C the opposite of CART3CMP and CART2CMP. The CMP2CART is used in bcktrn.
C The other is not used. ACES II gurus tend to use "computational" as
C some sort of ordering thing, but in here it only refers to the user's choice
C of Cartesian or Spherical basis sets.
C
C There is another convention called VMOL/ZMAT order. Obviously, for molecules
C with symmetry, coordinates of only the symmetry unique atoms are considered
C in the input. (The MOL files generated by JODA contain only the symmetry
C uniquie atoms) Now, both programs (VMOL and Seward) have internal logic to
C generate the non-unique atoms (coordinates) by applying the symmetry
C operations of the point group of the molecule. The order that these
C redundant atoms are generated by integral programs (VMOL or Seward) does
C not necessarily match with the the order that the user had choosen to enter
C them in the input file (ZMAT in ACES II, MOLCAS.INP in MOLCAS). This is the
C origin of the VMOL/ZMAT ordering. (For molecules with no symmetry this is a
C one to one map). Now, what I need to worry about is differences or
C similarities of the algorithms in Seward and VMOL that generate non-
C redundant centers. In VMOL these are done in readin.f (a scaled down
C version that I call from vmol2ja for seward runs to generate transformation
C matrices). In seward, these are done in rdctl.f, input.f, and their
C dependencies. The algorithms that are used for this purpose are very
C complicated and require a lot of patience to decipher (and, of course,
C a good knowledge of group theory. I hope Rod's new high-prized
C ACESQC/KDI/SBIR/GUI/... people can do this). A test program was written to
C test the assumption that seward and vmol generate symm. redundant atoms in
C the same order. Even though, I am not willing to bet my life on this,
C they indeed come out in the same order. One should not be suprised by this,
C since Roland had access to VMOL (MOLECULE in older days). His ideas about
C symmetry adaptation must have come from our old glorious friend. He has
C certainly done it a more modern way, but in my opinion the idea is the same.
C As I stated earlier, the algorithm is somewhat complicated and abstract,
C and their is no better way of doing than what it is in VMOL. 
C
C   MAP2ZMAT : This record contains the actual mapping from ZMAT to VMOL/SEWARD
C              (This is used in joda, symcor, libr, props, vscf)
C
C   CMP2ZMAT : This will transform VMOL/SEWARD ordered symm. adap. matrices to
C              ZMAT-ordered non sym. adapted Cartesian matrices. Here I am
C              purposefully avoiding the phrase "computational"
C              order. There is nothing computationally ordered
C              as far as I am concerned - only VMOL/SEWARD or ZMAT ordered
C              (correllated calculations are run in a totally different
C              order and are of no concern here). However, there is something
C              called computational basis (spherical or Cartesian depending
C              on user's choice). This transformation will work 
C              any computational basis. The ZMAT2CMP is the opposite
C              and these matrices are used in vscf, vscf_ks, intgrt, vksdint,
C              mrcc.
C
C The use of ZMAT2CMP and its inverse:
C    ZM=(ZMAT Order, non Symm. adapted)
C    VM=(VMOL or SEWARD Order, Symm. adapted)
C
C   CMP2ZMAT(ZM,VM)*MAT(VM,VM)*Trans(CMP2ZMAT) -> MAT(ZM,ZM)-A
C
C   ZMAT2CMP(VM,ZM)*MAT(ZM,ZM)*Trans(ZMAT2CMP) -> MAT(VM,VM)-B
C
C Note that A is in Cartesian non.-sym. adapted basis regardless of the
C computational basis of MAT(VM,VM) (spherical or cartesian), and B is in
C Computational sym. adapted basis regardless of the origin of MAT(ZM,ZM)
C (spherical or Cartesian).
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Common blocks
C     memman.inc  ICore   Dynamic core allocation
C     memini.inc  IStart  Index in ICore to begin allocation
C     mach.inc    MachSp  Machine-specific parameters
C     flags.inc   Flags   CRAPS & JODA control flags
C
C $Id: v2ja.F,v 1.6 2008/11/19 15:48:33 perera Exp $
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C     Include 'flags.inc'
C $Id: v2ja.F,v 1.6 2008/11/19 15:48:33 perera Exp $
      Integer IFlags(100), JFlags(16), Icore(Icrsiz)
      Common /Flags/ IFlags
      Common /Flags3/ JFlags
C
C Parameters
C
C     Include 'stdio.inc'
C $Id: v2ja.F,v 1.6 2008/11/19 15:48:33 perera Exp $
C     LuIn    Standard input
C     LuOut   Standard output
C     LuErr   Standard error output
C
      Integer LUIn, LuOut, LuErr
      Parameter (LuIn = 5, LuOut = 6, LuErr = 0)
C
C     LuInt   VMol or ``New VMol`` integral file
C     LuJobA  CRAPS JOBARC unit number
C
      Integer LuInt, LuJobA
      Parameter (LuInt = 1, LuJobA = 20)
C
C     FNVMol  VMol integral file name
C     FNNewV  New VMol integral file name
C     FNJobA  CRAPS JOBARC file name
C
      Character*(*) FNVMol, FNNewV, FNJobA
      Parameter (FNVMol = 'INT     ', FNNewV = 'IIII    ',
     &           FNJobA = 'JOBARC  ')
#include "maxbasfn.par"
C
C     ITNone  No integrals found
C     ITVMol  Integrals are from VMol
C     ITNewV  Integrals are from New VMol
C
      Integer ITNone, ITVmol, ITNewV, ITSew
      Parameter (ITNone = 0, ITVMol = 1, ITNewV = 2, ITSew = 3)
C
C     MxIrr   Maximum number of symmetry irreps possible
C
      Integer MxIrr
      Parameter (MxIrr = 8)
C
C Local Variables
C
C     IsTher  Used for checking existance of files
C     IntTyp  Typ of integrals being read
C     RepNuc  Nuclear repulsion energy
C     NBF     Nr. of (SO) basis functions
C     NBFIrr  Nr. (SO) basis fns. per irrep
C     NIrr    Nr. of irreps
C     NAOBF   Nr. of AO basis functions
C     Title   Title from integral file
C     i       General use
C     j       General use
C     LJunk   Logical for silly uses
C     MemPtr  Pointer array for memory allocation
C     MxMPtr  Max. number of MemPtr spaces to be used
C     MemUse  Max. memory used
C     NAtoms  Number of atoms
C     Det     Determinant from MINV
C     NOrbits Number of orbits (sym. unique atoms)
C     S       Scaling factor for computing inverse of CMP2ZMAT
C     NFAOIrr Number of functions per irrep in the full AO basis
C     SEWARD  If it is true, SEWARD integrals are being used.
C     ACESIII  If it is true, ACESIII is used.
C
      Double precision RepNuc, Det
      Integer IntTyp, NBF, NBFIrr(MxIrr), NIrr, i, MemPtr(15), MemUse,
     $   NAtoms, MxMPtr, NAOBF, NOrbits, LScr, IErr, IOrdGp,
     $   NFAOIrr(MxIrr), j, IACESIII
      Logical IsTher, LJunk, SEWARD, ACESIII
      Character*8 Title(24)
      CHARACTER*4 CENLBL(maxbasfn),ANGLBL(maxbasfn)
      Character*80 FNAME
C
C     **********************************************
C     * Setup to use CRAPS I/O and allocate memory *
C     **********************************************
C     By waiting until after we've checked for the existance of an
C     integral file, we may avoid wasting the effort of allocating
C     memory, etc.
C
CSSS      Call CrapsI (ICore(1), LJunk, 0)
C
      SEWARD = IFLAGS(56).EQ.4
      CALL GETREC(-1, "JOBARC", "AIIIJARC", 1, IACESIII)
      IF (IACESIII .GT. 0) ACESIII = .TRUE.
C
C     *************************
C     * Find an integral file *
C     *************************
C
      IntTyp = ITNone
C
      IF (SEWARD) Then
C
         Call Gfname('ONE_INT ',Fname, Ilength)
         Inquire (File=Fname(1:ILENGTH), Exist=IsTher)
C
         If (IsTher) then
            Open( LuInt, File=Fname(1:Ilength), Form='UNFORMATTED',
     $            Access='SEQUENTIAL')
            IntTyp = ITSew
         Endif
C
      Else If (ACESIII) Then
            IntTyp = ITNewV
      Else 
C
         Call Gfname(FNNewV,Fname,Ilength)
         Inquire (File=Fname(1:ILENGTH), Exist=IsTher)
C
         If (IsTher) then
            Open( LuInt, File=Fname(1:Ilength), Form='UNFORMATTED',
     $           Access='SEQUENTIAL')
            IntTyp = ITNewV
         Else
            Call Gfname(FNVMol,Fname,Ilength)
            Inquire (File=Fname(1:Ilength), Exist=IsTher)
            If (IsTher) then
               Open( LuInt, File=Fname(1:Ilength), Form='UNFORMATTED',
     $              Access='SEQUENTIAL')
               IntTyp = ITVMol
            EndIf
         EndIf
C
      EndIf
C
C     Make sure we found an integral file
C
      If ( IntTyp .eq. ITNone ) then
         Write (LuErr, 9500)
         Call ErrEx
      EndIf
 9500 Format (1X, '@V2JA-F, Unable to find integral file.')
C
C     ****************************
C     * Get the easy stuff first *
C     ****************************
C     Nuclear repulsion energy and nr. basis fns. per irrep
C
      If (.NOT. ACESIII) Then
         Rewind (LuInt)
         Read (LuInt) Title, NIrr, (NBFIrr(i), i=1, NIrr), RepNuc
C
C     Compute total number of (SO) basis fns
C
         NBF = 0
         Do 100 i = 1, NIrr
            NBF = NBF + NBFIrr(i)
 100     Continue
      Else
         Nirr = 1
         Call Getrec(20, "JOBARC", "NBASTOT ", 1, NBF)
         NBFIrr(1) = NBF
      Endif 
     
C     And the total number of AO basis functions
C
C     In principle, now that the CSYMTRAN record is written, we can
C     obtain NAOBF from it (below).  Maybe we we should eventually
C     get rid of AOSize.
C
      IF (.NOT. (SEWARD .OR. ACESIII)) THEN
         Call AOSize(LuInt, IntTyp .eq. ITNewV, NBF, NAOBF)
      ENDIF
C
C When SEWARD is used the SYMTRANS and CSYMTRANS records are
C not get written to the integral (IIII) file. Let's create
C them here to write them to the JOBARC file. Ajith 06/2000
C     
      
      If (SEWARD) Then
         Call GenTranMat(Icore, Icrsiz)
      endif
C
C     If the computational basis has dropped functions or it is spherical
C     harmmonic instead of cartesian, find out how many functions are in
C     each irrep of the full AO basis
C
      IF( .NOT. (SEWARD .OR. ACESIII))THEN
        If (IFlags(62) .eq. 1) then
           Call Locate(LuInt, 'CSYMTRAN')
           Read (LuInt) NAOBF, J, J, J, (NFAOIrr(i), i = 1, NIrr)
        EndIf
      ENDIF
C
C
C     ***************************************
C     * Write those records we know already *
C     ***************************************
C
      Call PutRec (LuJobA, FNJobA, 'NUCREP  ', IIntFP, RepNuc)
      Call PutRec (LuJobA, FNJobA, 'NBASTOT ',      1, NBF)
      Call PutRec (LuJobA, FNJobA, 'NIRREP  ',      1, NIrr)
      Call PutRec (LuJobA, FNJobA, 'NUMBASIR',   NIrr, NBFIrr)
cSSS      write(6,*) '  @V2JA-I, NAOBF is ',NAOBF
      IF(SEWARD .OR. ACESIII)THEN
cSSS        write(6,*) '  @V2JA-I, NAF is ',NBF
cSSS        write(6,*) '  @V2JA-I, NAOBF is ',NAOBF
        Call GETREC (LuJobA, FNJobA, 'NAOBASFN',      1, NAOBF)
cSSS        write(6,*) '  @V2JA-I, NAOBF is ',NAOBF
        Call GETREC (LuJobA, FNJobA, 'FAOBASIR',   NIrr, NFAOIrr)
cSSS        write(6,*) '  @V2JA-I, NFAOIrr ',NFAOIrr
      ELSE
        Call PutRec (LuJobA, FNJobA, 'NAOBASFN',      1, NAOBF)
        Call PutRec (LuJobA, FNJobA, 'FAOBASIR',   NIrr, NFAOIrr)
      ENDIF
cSSS      write(6,*) '  @V2JA-I, NAOBF is ',NAOBF
C
C     *******************************
C     * Allocate memory for GetAOSO *
C     *******************************
C     Need number of atoms for allocation; if its not there, abort
C
      I = 1
      Call GetRec (-LuJobA, FNJobA, 'NATOMS  ', I, NAtoms)
      I = 1
      Call GetRec( -1, 'JOBARC', 'COMPORDR', I, IORDGP)
C
C     Empty the pointer array and load it up with the lengths required
C
      MxMPtr = 15
      Call IZero(MemPtr, MxMPtr)
C
C     AO -> SO transformation matrix
C
      MemPtr( 1) = NAOBF * NAOBF * IIntFP
C
C     Atom and angular momentum labels; assume each integer can hold
C     a character*4.
C
c YAU - CenLbl() and AngLbl() replace MemPtr(2) and MemPtr(3)
c      MemPtr( 2) = NBF
c      MemPtr( 3) = NBF
      MemPtr( 2) = 0
      MemPtr( 3) = 0
C
C     NAOOrb, NSOOrb, IAngBF, NAOWrk, IAOCen
C
      MemPtr( 4) = NBF
      MemPtr( 5) = NAtoms
      MemPtr( 6) = NAOBF
      MemPtr( 7) = NAOBF
      MemPtr( 8) = NAOBF
C
C     NCenW1, NCenW2
C
      MemPtr( 9) = NAtoms
      MemPtr(10) = NAtoms
C
C     Coord1, Coord2
C
      MemPtr(11) = 3 * NAtoms * IIntFp
      MemPtr(12) = 3 * (IOrdGp-1) * IIntFP
C
C     CenMap
C
      MemPtr(13) = NAtoms
C
C     For debugging
C
      LScr = Max( NAOBF, NBF)
      MemPtr(14) = IIntFP * LScr * LScr
      MemPtr(15) = IIntFP * LScr * LScr
C
C     Turn the lengths into indices, start allocating a location 
C     returned by the dynamic memory allocator
C
      j = MemPtr(1) + Mod( MemPtr(1), 2)
      MemUse = j
      I0 = 1
      MemPtr(1) = I0
      Do 200 I = 2, MxMPtr
         iTmp = MemPtr(i) + Mod( MemPtr(i), 2)
         MemUse = MemUse + iTmp
         MemPtr(i) = MemPtr(i-1) + j
         j = iTmp
 200  Continue
C
C     Check if we've managed to overrun memory (unlikely, as our needs
C     are small, but still...
C
      If ( MemUse .gt. iCrSiz ) then
         call aces_fin
         Call InsMem( 'V2JA', MemUse, iCrSiz )
      EndIf
C
C     ****************
C     * Do some work *
C     ****************
C     Writes AO2SO, NBASATOM, CENTERBF, and ANGMOMBF
C
      Call Driver (LuInt, IntTyp .eq. ITNewV, NAOBF, NBF, NAtoms,
     $   NOrbits, IOrdGp,
     $   ICore(MemPtr(1)),  CenLbl,            AngLbl,
     $   ICore(MemPtr(4)),  ICore(MemPtr(5)),  ICore(MemPtr(6)),
     $   ICore(MemPtr(7)),  ICore(MemPtr(8)),  ICore(MemPtr(9)),
     $   ICore(MemPtr(10)), ICore(MemPtr(11)), ICore(MemPtr(12)),
     $   ICore(MemPtr(13)), LScr, ICore(MemPtr(14)), ICore(MemPtr(15)))
C
C     *****************************
C     * Allocate memory for Setup *
C     *****************************
C
C     Empty the pointer array and load it up with the lengths required
C
      MxMPtr = 12
      Call IZero(MemPtr, MxMPtr)
C
C     EVEC: Holds eigenvectors to be reordered (actually SOAO matrix)
C
      MemPtr( 1) = NAOBF * NAOBF * IIntFP
C
C     SCR: Scratch array for molecular coordinates & eigenvector
C
      LScr = Max (NAOBF, NBF, 3, NAtoms)
      MemPtr( 2) = LScr * LScr * IIntFP
C
C     NBASCN:
C
      MemPtr( 3) = NAtoms
C
C     COORD: Molecular coordinates
C
      MemPtr( 4) = NAtoms * 3 * IIntFP
C
C     IMEMCP, IPOPCMP, IANGCMP, IMAP, IBASOFF, IDUMMY
C
      MemPtr( 5) = NAOBF
      MemPtr( 6) = NAtoms
      MemPtr( 7) = NAOBF
      MemPtr( 8) = NAtoms
      MemPtr( 9) = 2 * NAtoms
      MemPtr(10) = NAOBF
C     For testing purposes
      MemPtr(11) = NAOBF
      MemPtr(12) = NAOBF
C
C     Turn the lengths into offsets, start allocating a location
C     returned by CrapsI
C
      Do 250 I = 2, MxMPtr
C
C        Make sure an even number of words is allocated for those
C        machines which require arrays to start on doubleword bounds
C
         MemPtr(i) = MemPtr(i) + Mod( MemPtr(i), 2)
C
C        Turn the specified length into an offset; but the offsets
C        are shifted one element below the desired variable.
C
         MemPtr(i) = MemPtr(i) + MemPtr(i-1)
 250  Continue
C
C     Correct the pointers being off by one position
C
C
      MemUse = MemPtr(MxMPtr) - I0
      Do 260 i = MxMPtr, 2, -1
         MemPtr(i) = MemPtr(i-1)
 260  Continue
      MemPtr(1) = I0
C
C     Check if we've managed to overrun memory (unlikely, as out needs
C     are small, but still...
C
      If ( MemUse .gt. iCrSiz ) then
         call aces_fin
         Call InsMem( 'V2JA', MemUse, iCrSiz )
      EndIf
C
C     *********************
C     * Do some more work *
C     *********************
C     Work on reordering the FULLSOAO matrix, already in ICore(Memptr(1))
C
C     Writes MAP2ZMAT, ANMOMBF0, CNTERBF0
C
C     Since we're sending in the full SOAO transformation, tell it
C     about the full (NAOBF instead of NBF) row dimension.
C

      Call MkC2Z1 (ICore(MemPtr(1)), ICore(MemPtr(2)), LScr,
     $   ICore(MemPtr(3)),
     $   ICore(MemPtr(4)), ICore(MemPtr(5)), ICore(MemPtr(6)),
     $   ICore(MemPtr(7)), ICore(MemPtr(8)),
     $   ICore(MemPtr(10)), NAOBF, NAtoms, ICore(MemPtr(11)), NAOBF )

C
C     If the basis is generally contracted, we need to rearrange the
C     p, d, etc. orbitals.
C
      I = 16
      Call GetRec (LuJobA, FNJobA, 'JODAFLAG', I, JFlags)
C
      If ( JFlags(14) .eq. 1) then
C
C        Get the l-angular momentum types
C
         Call GetRec(-LuJobA, FNJobA, 'ANMOMBF0', NAOBF,
     $      ICore(MemPtr(7)))
C
C        Do the rearrangement
C
         Call PDSwap (ICore(MemPtr(1)), ICore(MemPtr(7)),
     $      ICore(MemPtr(2)), NAOBF, NAOBF)
      EndIf
C
C     Want to write out the NAOBF x NBF submatrix of the full 
C     transformation.  Since that simply leavees off some columns
C     we don't have to do anything special.
C
      Call PutRec( LuJobA, FNJobA, 'CMP2ZMAT', IIntFP*NAOBF*NBF,
     $   ICore(MemPtr(1)) )
C
C     Find the inverse of this beastie
C
CSSS      call output(ICore(MemPtr(2)), 1, NAOBF, 1, NAOBF, NAOBF, NAOBF, 1)
      Call SCOPY(NAOBF*NAOBF, ICore(MemPtr(1)), 1, ICore(MemPtr(2)), 1)
c OLD
c      Call MInv(ICore(MemPtr(2)), NAOBF, NAOBF, ICore(MemPtr(1)), Det,
c     $   1.0d-8, 0, 1)
c NEW
      Call DGETRF(NAOBF,NAOBF,ICore(MemPtr(2)),NAOBF,ICore(MemPtr(1)),
     &            iErr)
      iErr = 0
      if (iErr.ne.0) then
         if (iErr.gt.0) then
            print *, '@V2JA: factorized a singular matrix'
            print *, '       Change the basis set.'
         else
            print *, '@V2JA: dgetrf argument ',-iErr,' is illegal'
         end if
         call errex
      end if
      iPad = iand(NAOBF,ishft(iand(IINTFP,2),-1))

      Call DGETRI(NAOBF,ICore(MemPtr(2)),NAOBF,ICore(MemPtr(1)),
     &            ICore(MemPtr(1)+NAOBF),NAOBF*(NAOBF-1),iErr)
      if (iErr.ne.0) then
         if (iErr.gt.0) then
            print *, '@V2JA: inverted a singular matrix'
            print *, '       Change the basis set.'
         else
            print *, '@V2JA: dgetri argument ',-iErr,' is illegal'
         end if
         call errex
      end if
c END
      IErr = 0
cSSS      Call XGEMM('N', 'N', NAOBF, NAOBF, NAOBF, 1.0d0, ICore(MemPtr(2)),
cSSS     $   NAOBF, ICore(MemPtr(1)), NAOBF, 0.0d0, ICore(MemPtr(3)), NAOBF)
cSSS      Call DiagCk(NAOBF, NAOBF, 1.0d-8, ICore(MemPtr(3)), NAOBF, LScr)
cSSS      Write (6, *) 'Off-diagonals = ', LScr
C
C     Write out the NBF x NAOBF submatrix of it
C
      Call BlkCpy2(ICore(MemPtr(2)), NAOBF, NAOBF,
     $   ICore(MemPtr(1)), NBF, NAOBF, 1, 1)
cSSS      Call XGEMM('N', 'N', NAOBF, NAOBF, NAOBF, 1.0d0, ICore(MemPtr(3)),
cSSS     $   NAOBF, ICore(MemPtr(1)), NAOBF, 0.0d0, ICore(MemPtr(2)), NAOBF)
cSSS      Call DiagCk(NAOBF, NAOBF, 1.0d-8, ICore(MemPtr(2)), NAOBF, LScr)
cSSS      Write (6, *) 'Off-diagonals = ', LScr
      Call PutRec( LuJobA, FNJobA, 'ZMAT2CMP', IIntFP*NBF*NAOBF,
     $   ICore(MemPtr(1)) )
C
C     ************
C     * All done *
C     ************
C
      Return
      End
