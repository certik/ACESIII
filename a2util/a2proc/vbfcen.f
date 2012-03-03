C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine VBFCen(NOrb, CenOrb, NBFOrb, BF2Cen)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: vbfcen.f,v 1.1.1.1 2003/04/02 19:21:47 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     VBFCen -- generate basis function to center mapping for VMol
C
C SYNOPSIS
      Integer NOrb, CenOrb(NOrb), NBFOrb(NOrb), BF2Cen(*)
C
C ARGUMENTS
C     NOrb    Number of orbits, aka unique centers (input)
C     CenOrb  Number of center in each orbit (input)
C     NBFOrb  Number of basis fns in each orbit (input)
C     BF2Cen  Mapping of basis function into center number (output)
C
C DESCRIPTION
C     Uses an understanding of the innards of VMol to create a mapping
C     of the AO basis set to the corresponding center.
C
C NOTES
C     NBFOrb is the number of basis functions for the whole orbit,
C     not the number of functions on a given center in the orbit.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     BF      Counter for basis functions
C     Cen     Cummulative number of centers before the current orbit
C     I       Loop over functions in orbit
C
      Integer BF, Cen, I, IOrb, j
C
C     Initialize counters
C
      BF = 1
      Cen = 0
C
C     Loop over the orbits, aka unique centers
C
      Do 100 IOrb = 1, NOrb
C
C        Place one copy of a basis function on each image equivalent
C        center before moving on to the next basis function.
C
C        Loop over the number of basis fns on each equiv. center.
C
         Do 110 I = 1, ( NBFOrb(IOrb) / CenOrb(IOrb) )
C
C           Loop through the equiv. centers in order to give each
C           a copy of this basis fn.
C
            Do 120 J = 1, CenOrb(IOrb)
               BF2Cen(BF) = Cen + J
               BF = BF + 1
 120        Continue
 110     Continue
C
C        Move Cen up to the beginning of the next orbit
C
         Cen = Cen + CenOrb(IOrb)
 100  Continue
      Return
      End
