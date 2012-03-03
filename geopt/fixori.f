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
      Subroutine FIXORI (NAtms, Q, At1, At2, OldOri)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Integer NAtms, At1, At2, OldOri
      Dimension Q(3*NAtms),ORIENT(3,3),ORIENT2(3,3)
C
      Integer NewOri
      Double Precision NewVec(3), Tmp
C     Define the "Standard I/O" units to use throughout the code
      Parameter (LuIn  = 5)
      Parameter (LuOut = 6)
      Parameter (LuErr = 6)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C     Its quite possible that there won't be anything to check...
C
      If (At1 .eq. 0 .AND. At2 .eq. 0 .AND. OldOri .eq. 0) Return
C
C     Get the vector from At1 to At2 and find out how it is directed
C
      Call VAdd(NewVec, Q(3*(At2-1)+1), Q(3*(At1-1)+1),3,-1.D0)
      Call Normal(NewVec,3)
      NewOri = NInt(NewVec(1)) + 2*NInt(NewVec(2))
     $   + 3*NInt(NewVec(3))
C
C     Check to see if this orientation is the same as last time.  If
C     not, interchange x and y axes, keeping the coordinate frame
C     right-handed.
C
C     Note:  Handedness of frame may change eariler, so atoms may now
C     be connected by the *opposite* vector.  This is okay.
C
      If (Abs(NewOri) .ne. Abs(OldOri)) then
         ININE=9
         CALL DGETREC(20,'JOBARC','ORIENT2 ',ININE,ORIENT2)
         CALL DGETREC(20,'JOBARC','ORIENTMT',ININE,ORIENT)
         Write (LuOut, 9000)
 9000    Format ('@FIXORI-I, Rotating about z-axis to ',
     $      'match previous orientation.')
         Do 100 i = 1, NAtms
            Tmp = Q(3*(i-1)+1)
            Q(3*(i-1)+1) = Q(3*(i-1)+2)
            Q(3*(i-1)+2) = -Tmp
 100     Continue
         TMP=ORIENT(1,1)
         ORIENT(1,1)=ORIENT(2,2)
         ORIENT(2,2)=TMP
         TMP=ORIENT(1,2)
         ORIENT(1,2)=ORIENT(2,1)
         ORIENT(2,1)=TMP
         TMP=ORIENT2(1,1)
         ORIENT2(1,1)=ORIENT2(2,2)
         ORIENT2(2,2)=TMP
         TMP=ORIENT2(1,2)
         ORIENT2(1,2)=ORIENT2(2,1)
         ORIENT2(2,1)=TMP
         CALL DPUTREC(20,'JOBARC','ORIENT2 ',ININE,ORIENT2)
         CALL DPUTREC(20,'JOBARC','ORIENTMT',ININE,ORIENT)
         CALL DPUTREC(20,'JOBARC','COORD   ',NATMS*3,Q)
      EndIf
      Return
      End
