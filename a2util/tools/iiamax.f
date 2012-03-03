
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     anthony yau, modified directly from idamax

      integer function iiamax(n,iarray,incx)
      integer n,iarray(1),incx
      integer i,ix,imax,itmp

      if (n.le.1) then
         if (n.eq.1) then
            iiamax=1
         else
            iiamax=0
         end if
         return
      end if

      itmp=1
      imax=abs(iarray(1))
      if (incx.eq.1) then
         do i=2,n
            if (abs(iarray(i)).gt.imax) then
               itmp=i
               imax=abs(iarray(i))
            end if
         end do
      else
         ix=1
         do i=2,n
            ix=ix+incx
            if (abs(iarray(ix)).gt.imax) then
               itmp=i
               imax=abs(iarray(ix))
            end if
         end do
      end if
      iiamax=itmp

      return
      end

