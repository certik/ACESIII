
C     mkc2z1 -- First step in producing the transformation from the
C     computational basis to the "nice" ZMAT-ordered basis
C
C DESCRIPTION
C     Transforms EVEC from the VMol AO basis (not SO basis!) to the
C     "ZMAT" order preferred for symmetry analysis.  Primarily, this
C     means changing the order of the atoms.  Also produces reordered
C     versions of CENTERBF and ANGMOMBF on the JOBARC.
C
C NOTES
C     Some blocks of comments in the code have been marked with "C--"
C     before and after them.  These were the original comments put in
C     by the author.  Other comments have been added in an attempt to
C     understand what the heck is going on in this routine!

      SUBROUTINE MkC2Z1(EVEC,SCR,LSCR,NBASCN,COORD,IMEMCMP,IPOPCMP,
     &                  IANGCMP,IMAP,IDummy,NBAS,NATOMS,ICheck,NAO)

      Integer NAO, LScr, NBas, NAtoms
      Double precision EVEC(NAO,NBAS), SCR(LSCR,LSCR), Coord(3,NAtoms)
      Integer NBASCN(NATOMS),IMEMCMP(NAO),IPOPCMP(NATOMS),
     $   IMAP(NATOMS),IANGCMP(NAO), IDUMMY(NAO)
      Integer ICheck(NAO)
C     LSCR should be Max(NAO, NBas, 3, NAtoms)

      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      Integer I, IOne, IThr, IOrdGp, NOrbit
      CHARACTER*4 CMPPGP

#ifdef _TRACE_VMOL2JA
      write(*,*) '@VMOL2JA: Entered ', __FILE__
#endif /* _TRACE_VMOL2JA */

      IONE=1
      ITHR=3*IINTFP
C
C--
C PICK UP SOME INFORMATION FROM THE JOBARC FILE
C--
C
      CALL GETREC(-1,'JOBARC','COMPNORB',IONE,NORBIT)
      CALL GetREC(20,'JOBARC','MAP2ZMAT',NATOMS,IMAP)
      Call GetRec(20,'JOBARC','CENTERBF',NAO,IDummy)
C
C     Renumber the centers from VMol to ZMAT arrangement
C
      Do 10 I = 1, NAO
         IDummy(i) = IMap( IDummy(i) )
 10   Continue
C
C     John's preferred ordering
C
      Call JFSOrd(NAtoms, NAO, IDummy, ICheck)

C     Reorder CENTERBF, already in core.
      Call ISctr(NAO,IDummy,ICheck,IAngCmp)
#ifdef _DEBUG_V2JA
      write(*,*) '@MKC2Z1: Recording CNTERBF0 as '
      write(*,*) '         ', (IAngCmp(i), i=1,NAO)
#endif /* _DEBUG_V2JA */
      Call PutRec(20,'JOBARC','CNTERBF0',NAO,IAngCmp)

C     Reorder ANGMOMBF
      Call GetRec(20,'JOBARC','ANGMOMBF',NAO,IDummy)
      Call ISctr(NAO,IDummy,ICheck,IAngCmp)
#ifdef _DEBUG_V2JA
      write(*,*) '@MKC2Z1: Recording ANMOMBF0 as '
      write(*,*) '         ', (IAngCmp(i), i=1,NAO)
#endif /* _DEBUG_V2JA */
      Call PutRec(20,'JOBARC','ANMOMBF0',NAO,IAngCmp)

cSSS      Call GetRec( 20, 'JOBARC', 'CNTERBF0', NAO, IDummy)
cSSS      Call NewPD (NAO, IDummy, IAngCmp, IMemCmp, I)
cSSS      Write (6, 9000) 'NEWPD', (i, IMemCmp(i), i = 1, NAO)
C
C     Reorder the transformation; written back out by the caller.
C
      Do 100 I = 1, NAO
         Call xCOPY(NBas, Evec(i, 1), NAO, Scr( ICheck(i), 1), NAO)
 100  Continue
      Call xCOPY(NAO*NBas, Scr, 1, Evec, 1)
C
 9000 Format(1X,'new style ', A,':'/(2I5))

#ifdef _TRACE_VMOL2JA
      write(*,*) '@VMOL2JA: Leaving ', __FILE__
#endif /* _TRACE_VMOL2JA */

      RETURN 
      END
