C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      SUBROUTINE AOSize(LUN, NewVM, NSO, NrAO)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: aosize.f,v 1.5 2008/11/19 15:48:33 perera Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     aosize -- determine size of the AO basis for VMol
C
C SYNOPSIS
      Integer LUN, NSO, NrAO
      Logical NewVM
C
C ARGUMENTS
C     LUN     Unit number of integral file (input)
C     NSO     Nr fo symmetry-adapted basis functions (input)
C     NewVM   New VMol or old VMol (input)
C     NrAO    Largest AO index referenced in AO->SO xform. (output)
C
C ROUTINES REQUIRED
C     IAMAX   Largest absolute element of a vector
C     Locate  Find a Molecule-style record label 
C
      External IAMAX, Locate
      Integer IAMAX
C
C DESCRIPTION
C     Determines the size of the AO basis by finding the largest function
C     index referenced by the AO->SO transformation.  Anything higher
C     than that can't possibly be relevant... we hope!
C
C BUGS
C     Limited to 80 AOs contributing to each SO.
C     Molecule record label is fixed as SYMTRANS.
C     Variable IIntFP in common block MACHSP must be 1 or 2.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Include 'mach.inc'
C $Id: aosize.f,v 1.5 2008/11/19 15:48:33 perera Exp $
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C LOCAL VARIABLES
C     i       General counter
C     ii      General counter
C     j       General counter
C     JTran   Number of AOs in SO (valid elements in ITran/CTran)
C     ITran   AO components of an SO
C     Junk    To skip over center and ang. mom. labels neatly.
C     SOAO    Scratch array to read coefficients into
C             
      Character*4 Junk
      Integer I, ii, j, JTran, ITran(80)
      Double Precision SOAO
C
C     Find the symmetry transformation section
      Call Locate( LUN, 'SYMTRANS')
C
C     Record format is slightly different for New_Mol & VMol
C        New_Mol starts with a record holding the number of basis fns
C           and is otherwise like VMol.
C        VMol writes each orbital as a record
C
C     NOTE:  The number of 4 char. pieces making up the two labels
C     should be 4 if IIntFP=1 and 2 if IIntFP=2.
C
C
      NrAO = 0
      If (.NOT. NewVM) then
         Do 100 I = 1, NSO
            Read (LUN) J, (Junk, ii = 1, 2*(3-IIntFP)),  JTran,
     $         (ITran(ii), SOAO, ii = 1, JTran)
               NrAO = Max( NrAO, IAMax(JTran, ITran, 1) )
 100     Continue
      Else
         Read (LUN)
         Do 200 I = 1, NSO
            Read (LUN) JTran, (ITran(ii), SOAO,
     $         ii = 1, JTran)
               NrAO = Max( NrAO, IAMax(JTran, ITran, 1) )
 200     Continue
      EndIf
C
      Return
      End
