
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine RdSOAO (LUN, RecNam, NewVM, NAO, NSO, SOAO, LD)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: rdsoao.f,v 1.4 2008/11/19 15:48:33 perera Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rdsoao -- read SO --> AO transformation from MOLECULE file
C
C SYNOPSIS
      Integer LUN, NSO, NAO, LD
      Character*(*) RecNam
      Double precision SOAO(LD, NSO)
      Logical NewVM
C
C ARGUMENTS
C     LUN     Unit for integral file (input)
C     RecNam  Record name for transformaiton (input)
C     NewVM   New VMol or old VMol (input)
C     NAO     Row dimension SOAO (nr of AOs) (input)
C     NSO     Column dimension of SOAO (nr of SOs) (input)
C     SOAO    The SO -> AO transformation (output)
C     LD      Leading dimension of SOAO (input)
C
C ROUTINES REQUIRED
C     Zero    Zero out a vector
C     Locate  Find a give record label in a Molecule-style int. file
C
C DESCRIPTION
C     Reads SO -> AO transform. matrix from a MOLECULE integral file.
C     In the output matrix, SOs label the columns and AOs the rows.
C     
C     The usual record to read is labeled SYMTRANS.
C
C     The matrix is not in general normalized, so the inverse
C     transformation requires some playing with factors.  Moreover, 
C     the matrix may not be square if spherical harmonic SOs are
C     used.  In this case the matrix may still be considered square
C     with "phantom" columns representing the SOs which have been
C     dropped in the cartesian->spherical transformation.  This
C     expanded matrix would then have to act on a matrix with the SO
C     space augemented by rows of zeros corresponding to the dropped
C     functions.  Good luck!
C
C BUGS
C     Limited to 80 AOs contributing to each SO.
C     Variable IIntFP in common block MACHSP must be 1 or 2.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Include 'mach.inc'
C $Id: rdsoao.f,v 1.4 2008/11/19 15:48:33 perera Exp $
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C LOCAL VARIABLES
C     i       General counter
C     ii      General counter
C     j       General counter
C     JTran   Number of AOs in SO (valid elements in ITran/CTran)
C     ITran   AO components of an SO
C     Junk    To skip over center and ang. mom. labels neatly.
C             
      Character*4 Junk
      Integer I, ii, j, JTran, ITran(224)
C
      Call Zero( SOAO, LD*NSO)
C
C     ***************************
C     * SO -> AO transformation *
C     ***************************
C
C     Find the symmetry transformation section
      Call Locate( LUN, RecNam)
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
      If (.NOT. NewVM) then
         Do 100 I = 1, NSO
            Read (LUN) J, (Junk, ii = 1, 2*(3-IIntFP)),  JTran,
     $         (ITran(ii), SOAO(ITran(ii), j), ii = 1, JTran)
 100     Continue
      Else
         Read (LUN)
         Do 200 I = 1, NSO
            Read (LUN) JTran, (ITran(ii), SOAO(ITran(ii), i),
     $         ii = 1, JTran)
 200     Continue
      EndIf
C
      Return
      End
