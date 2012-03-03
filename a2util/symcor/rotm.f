
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

