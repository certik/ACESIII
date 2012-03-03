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

C     FORMS THREE DIMENSIONAL ROTATION MATRIX
C     (IX=0 IF ANG IN RADS,
C      IX=1 IF DEG,
C      ADD 10 TO IX TO GET TRANSPOSE)

      SUBROUTINE ROTM( IAXIS, ANG, IX, RT )
      implicit double precision (a-h,o-z)

      integer iaxis, ix
      dimension rt(3,3)

      integer iaxis1, iaxis2

      CALL ZERO(RT,9)
      IF (mod(IX,10).EQ.1) THEN
         TANG=ANG*DACOS(-1.D0)/180.D0
      ELSE
         TANG=ANG
      END IF

      if (ix/10.gt.0) then
         sign = -1.d0
      else
         sign = 1.0d0
      end if

c     next two axes
c     --  really mod( (iaxis-1)+1, 3 ) + 1 (clock, not mod, arithmetic)
      iaxis1 = mod(iaxis ,3)+1
      iaxis2 = mod(iaxis1,3)+1

      rt(iaxis,iaxis) = 1.0 d0
      RT(iaxis1,iaxis1) = COS(TANG)
      RT(iaxis2,iaxis2) = COS(TANG)

      rt(iaxis1,iaxis2) = sin(tang) * sign
      rt(iaxis2,iaxis1) = -rt(iaxis1,iaxis2)

      RETURN
      END

