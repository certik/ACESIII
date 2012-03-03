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
      Integer Function IOrGrp ( PtGrp )
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Purpose:
C     Return the order of the given point group
C
C Arguments:
C     PtGrp    String containing the Schoenflies symbol for the
C              group of interest. (input)
C              The order must begin with the second character of PtGrp.
C              The last non-blank character is taken to be the "modifier"
C              in the group symbol. It may be in upper- or lowercase.
C              A group of infinite order is denoted with an 'X' for the
C              order.
C
C Returns:
C     IOrGrp   Returns the order of the group.
C              In case of error, returns error codes a follows:
C              = -1  Unrecognizable point group
C
C Dependents:
C     AtoI     Converts string to integer
C     linblnk   Returns the index of the last non-blank character
C
C Bugs:
C     Since most machines can't actually handle a value of infinity,
C     infinite groups return >= 1,000,000 - which should be suitably
C     large to clue someone in.  Recall that Ih is order 120.
C
C     Doesn't check the modifier against the group - thus things like
C     C4d won't be flagged as an error.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Implicit Integer (A-Z)
      Character*(*) PtGrp
      Integer AtoI, linblnk
C
C     PGMod is the "modifier" lifted from the name
C
      Character*1 PGMod
C
C     Get the "modifier" from the end of the symbol
C
      i = linblnk(PtGrp)
      PGMod = PtGrp(i:i)
C
C     Figure out the order of the axis.  If there is no number 
C     there, AtoI returns 0
C
      OrdAx = AtoI( PtGrp(2:) )
C
C     Assume the symbol is bad.  If it isn't, it will be changed
C     as we go through everything.  Otherwise, this is the right
C     thing to return.
C     
      IOrGrp = -1
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Now check for the cubics & other groups without explicit orders
C     in the Schoenflies symbol.
C
      If (OrdAx .eq. 0) then
         If (PtGrp(1:1) .eq. 'T') then
            IOrGrp = 12
         ElseIf (PtGrp(1:1) .eq. 'O') then
            IOrGrp = 24
         ElseIf (PtGrp(1:1) .eq. 'I') then
            IOrGrp = 60
         ElseIf (PtGrp(1:1) .eq. 'K') then
            IOrGrp = 1000000
         ElseIf (PtGrp(1:2) .eq. 'CX' .OR. PtGrp(1:2) .eq. 'DX') then
            IOrGrp = 1000000
         ElseIf (PtGrp(1:1) .eq. 'C' .AND. (PGMod .eq. 's' .OR.
     $         PGMod .eq. 'S' .OR. PGMod .eq. 'i' .OR.
     $         PGMod .eq. 'I')) then
            IOrGrp = 2
         EndIf
      Else
C
C     For Cn and Sn, the order is n.  For Dn, it is 2n.
C
         If (PtGrp(1:1) .eq. 'C' .OR. PtGrp(1:1) .eq. 'S') then
            IOrGrp = OrdAx
         ElseIf (PtGrp(1:1) .eq. 'D') then
            IOrGrp = 2 * OrdAx
         EndIf
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Now handle the possible modifiers:
C     v, h, & d are all 2x the number for the base group.
C
      If (PGMod .eq. 'v' .OR. PGMod .eq. 'V' .OR. PGMod .eq. 'h'
     $   .OR. PGMod .eq. 'H' .OR. PGMod .eq. 'd' .OR. PGMod .eq. 'D')
     $   then
         IOrGrp = 2 * IOrGrp
      EndIf
C
C     Make sure we return the right error code.
C
      If (IOrGrp .lt. 0) IOrGrp = -1
      Return
      End
