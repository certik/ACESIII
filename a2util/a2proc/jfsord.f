C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine JFSOrd (NCen, NBF, BF2Cen, Map)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: jfsord.f,v 1.1.1.1 2003/04/02 19:21:47 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     jfsord -- Make permutation from original order to the "jfs" order
C
C SYNOPSIS
      Integer NCen, NBF
      Integer BF2Cen(NBF), Map(NBF)
C
C ARGUMENTS
C     NCen    Number of centers (input)
C     NBF     Number of basis functions (input)
C     Bf2Cen  Mapping of basis functons to centers (input)
C     Map     Mapping from BF2Cen order to "jfs" order (output)
C
C DESCRIPTION
C     The "jfs" order has all functions on a given center grouped
C     together and the thing in order by center number.  Care is taken
C     to insure that individual basis functions stay in the same
C     relative order.
C
C     The "jfs" order normally implies that the center numbers on which
C     this arrangement is based are the ZMAT center numbers. This
C     routine takes no explicit account of that, but will give the
C     desired result if the Bf2Cen array contains the ZMAT center nrs.
C
C     The result of this routine is Map, which should be used as:
C     Original(i) --> Reordered( Map(i) ).  This is conceptually a
C     "scatter" operation when applied to a vector.
C
C BUGS
C     Does not check for out-of-range center numbers in BF2Cen.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     IOff    Offset into Map
C     I       counter over centers
C     J       counter over basis
C
      Integer IOff, I, J
C
      IOff = 1
C
C     Look for each center in sequence
C
      Do 100 I = 1, NCen
C
C        Look for every occurrance of this center and use Map to collect 
C        then together into a contiguous block.
C
         Do 110 J = 1, NBF
            If (BF2Cen(j) .eq. i) then
               Map(j) = IOff
               IOff = IOff + 1
            EndIf
 110     Continue
 100  Continue
      Return
      End
