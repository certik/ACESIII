



C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine RdBLab (LUN, NewVM, NSO, Cen, AngMom)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: rdblab.F,v 1.2 2004/07/30 02:20:08 yau Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rdblab -- read basis function labels from Molecule-style integrals
C
C SYNOPSIS
      Integer LUN, NSO
      Character*4 Cen(NSO), AngMom(NSO)
      Logical NewVM
C
C ARGUMENTS
C     LUN     Unit for integral file (input)
C     NewVM   New VMol or old VMol (input)
C     NSO     Nr. of symmetry adapted basis functions (input)
C     Cen     SO function center labels (output)
C     AngMom  SO function angular momentum labels (output)
C
C ROUTINES REQUIRED
C     Locate  Find a give record label in a Molecule-style int. file
C
C DESCRIPTION
C     Reads the center an angular momentum text labels from a Molecule-
C     style integral file.
C
C     Center labels really identify the center selected to be the origin
C     of each orbit (set of equivalent centers).
C
C BUGS
C     Molecule record label is fixed as SYMTRANS for old VMol and
C     LABBASIS for new VMol.
C     Variable IIntFP in common block MACHSP must be 1 or 2.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Include 'mach.inc'
C $Id: rdblab.F,v 1.2 2004/07/30 02:20:08 yau Exp $
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
C LOCAL VARIABLES
C     i       General counter
C     j       General counter
C     Junk    To skip over center and ang. mom. labels neatly.
C             
      Character*4 Junk
      Integer I, j, IACESIII
      LOGICAL SEWARD, ACESIII
C
      SEWARD = IFLAGS(56) .EQ. 4
      CALL GETREC(-1, "JOBARC", "AIIIJARC", 1, IACESIII)
      IF (IACESIII .GT. 0) ACESIII = .TRUE.

      IF (SEWARD .OR. ACESIII) THEN
         OPEN(UNIT=1,FILE='LABBASIS',STATUS='OLD')
         DO I = 1,NSO
            READ(1,'(3X,A4,3X,A4)') CEN(I),ANGMOM(I)
         END DO
         CLOSE(UNIT=1,STATUS='KEEP')
         RETURN
      END IF
C
C     New and old VMol put the data in different places
C
      If (NewVM) then
         Call Locate( LUN, 'LABBASIS')
      Else
         Call Locate( LUN, 'SYMTRANS')
      EndIf
C
C     Also have to worry about padding due to word length considerations
C
C SG 8/26/98 Not anymore, since these are written as character*4
      Do 200 I = 1, NSO
         Read (LUN) J, Cen(j), AngMom(j)
 200  Continue
C
      Return
      End
