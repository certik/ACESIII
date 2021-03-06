C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Integer Function IAMAX( N, IX, IncX)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: iamax.f,v 1.1.1.1 2003/04/02 19:21:47 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     iamax -- maximum absolute value of an integer vector
C
C SYNOPSIS
      Integer N, IX(1), IncX
C
C ARGUMENTS
C     N       Number of elements to check (input)
C     IX      Integer vector (input)
C     IncX    Skip distance between elements of IX (input)
C
C DESCRIPTION
C     Returns the maximum absolute value of the vector.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     IPtr    Pointer to current element of IX
C     j       Counter over number of "actIXe" elements
C
      Integer IPtr, j
C
C     Initialize the max
C
      IAMax = 0
C
C     Get out if there is no work
C
      If (N .le. 0) Return
C
C     Handle increment = 1 differently
C
      If (IncX .eq. 1) then
         Do 100 IPtr = 1, N
            IAMax = Max(IAMax, Abs(IX(IPtr)))
 100     Continue
      Else
         IPtr = 1
         If (IncX .lt. 0) IPtr = (-N+1) * IncX + 1
         Do 200 j = 1, N
            IAMax = Max(IAMax, Abs(IX(IPtr)))
            IPtr = IPtr + IncX
 200     Continue
      EndIf
C
      Return
      End
