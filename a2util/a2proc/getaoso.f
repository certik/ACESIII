


C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      SUBROUTINE GETAOSO(NVMol, NAO, NBas, AOSO, ATMLBL, ANGLBL, IPopVC,
     $   IANGBF, IAOOrb, NCen, IOrbit, IORBPOP)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: getaoso.f,v 1.2 2005/08/03 14:39:08 yau Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     getaoso -- create the AO to SO and other transformations
C
C SYNOPSIS
      Logical NVMol
      Integer NAO, NBas, NCen, IOrbit
      CHARACTER*4 ATMLBL(NBAS),ANGLBL(NBAS)
      Double precision AOSO(NAO,NBAS)
      Integer IPopVC(NBas), IAngBF(NAO), IAOOrb(NAO), IOrbPOP(NCen)
C
      Integer Index
      Intrinsic Index
C
C DESCRIPTION
C     Reads a VMol integral file to obtain the necessary information
C     to create transformations or mapping functions for a number of
C     relations needed for symmetry analysis.  Results are written
C     to the CRAPS JOBARC file.
C
C     AO2SO     Transformation matrix from symmetry-adapted orbitals
C               to raw atomic orbitals .
C     NBASATOM  Number of SO basis functions per orbit.
C     NAOBFORB  Number of AO basis functions per orbit.
C     CENTERBF  Atomic center corresponding to each basis fn.
C     ^^^^^^^^  Moved to MkC2Z1
C     ANGMOMBF  l quantum number of each basis fn.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C
      Integer I, J, K
      Character*7 up_SPDFGHI, lo_spdfghi, XYZxyz
      up_SPDFGHI = 'SPDFGHI'
      lo_spdfghi = 'spdfghi'
      XYZxyz     = 'XYZxyzX'
C
      DO 119 I=1,NAO
         IANGBF(I)=-1
 119  CONTINUE
C
      DO 120 I=1,NBAS
        DO 121 J=1,NAO
           If (AOSO(J, I) .ne. 0.0d0) then
              IAngBF(J) = -1 + INDEX(up_SPDFGHI,AngLbl(i)(1:1))
     &                       + INDEX(lo_spdfghi,AngLbl(i)(1:1))
              IF (IAngBF(J).EQ.-1) THEN
                 IF (INDEX(XYZxyz,AngLbl(i)(1:1)).NE.0) IAngBF(J)=1
                 IF (INDEX(XYZxyz,AngLbl(i)(2:2)).NE.0) IAngBF(J)=2
                 IF (.NOT.NVMOL) THEN
                    IF (INDEX(XYZxyz,AngLbl(i)(3:3)).NE.0) IAngBF(J)=3
                    IF (INDEX(XYZxyz,AngLbl(i)(4:4)).NE.0) IAngBF(J)=4
                 END IF
              ENDIF
           EndIf
 121    CONTINUE
 120  CONTINUE
C
      CALL IZERO(IPOPVC,NBAS)
      CALL IZERO(IORBPOP, NCen)
      Call IZERO( IAOOrb, NAO)
C
C     Loop over basis functions
C
      IORBIT=0
      DO 130 I=1,NBAS
         IF(IPOPVC(I).NE.0) GOTO 130
         IORBIT=IORBIT+1
C
C        Compare each of the remaining atom labels with the one for the
C        outer loop.  If the labels match, this is another center for
c        the orbit.
C
         DO 131 J=I,NBAS
            IF(ATMLBL(I).EQ.ATMLBL(J)) THEN
               IPOPVC(J)=IORBIT
               IORBPOP(IORBIT)=IORBPOP(IORBIT)+1
C              
C              Mark the AOs that contribute to this SO as being in this
C              orbit too.
C              
               Do 200 K = 1, NAO
                  If (AOSO(K, J) .ne. 0.0d0) IAOOrb( K ) = IOrbit
 200           Continue
C              
            ENDIF
 131     CONTINUE
 130  CONTINUE
C
C     Now that we've assigned each of the AOs to an orbit, collect the
C     results into a more compact form.
C
      Call IZERO(IPOPVC, IOrbit)
      Do 220 I = 1, NAO
         IPopVC( IAOOrb(i) ) = IPopVC( IAOOrb(i) ) + 1
 220  Continue
C
      RETURN
      END
