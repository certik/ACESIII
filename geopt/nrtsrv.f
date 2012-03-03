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
      Integer Function NrTSRV ( PtGrp )
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Purpose:
C     Return the number of rotations and vibrations in the totally
C     symmetric irrep.
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
C Dependents:
C     linblnk   Returns the index of the last non-blank character
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Implicit Integer (A-Z)
      Character*(*) PtGrp
      Integer linblnk
C
C     PGMod is the "modifier" lifted from the name
C
      Character*1 PGMod
C
C     Get the "modifier" from the end of the symbol
C
      i = linblnk(PtGrp)
      PGMod = PtGrp(i:i)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     D groups and cubic groups
C
      If (PtGrp(1:1) .eq. 'D' .OR. PtGrp(1:1) .eq. 'T' .OR.
     $   PtGrp(1:1) .eq. 'O' .OR. PtGrp(1:1) .eq. 'I' .OR.
     $      PtGrp(1:1) .eq. 'K') then
         NrTSRV = 0
C
C     S groups
C
      ElseIf (PtGrp(1:1) .eq. 'S') then
         NrTSRV = 1
C
C     C groups: C1 => 6; Cs & Ci => 3; Cnh & Cnv => 1; Cn => 2
C
      ElseIf (PtGrp(1:1) .eq. 'C') then
         If (PtGrp(2:linblnk(PtGrp)) .eq. '1') then
            NrTSRV = 6
         ElseIf (PGMod .eq. 's' .OR. PGMod .eq. 'S' .OR. PGMod .eq. 'i'
     $      .OR. PGMod .eq. 'I') then
            NrTSRV = 3
         ElseIf (PGMod .eq. 'v' .OR. PGMod .eq. 'V' .OR. PGMod .eq. 'h'
     $      .OR. PGMod .eq. 'H') then
            NrTSRV = 1
         Else
            NrTSRV = 2
         EndIF
C
C     It isn't a group we recognize
C
      Else
         NrTSRV = -1
      EndIf
C
      Return
      End
