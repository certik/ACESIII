
c     finds the index of element having min. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     anthony yau, modified directly from idamin

      integer function iiamin(n,iarray,incx)
      integer n,iarray(*),incx
      integer i,ix,imin,itmp

      if (n.lt.2) then
         if (n.eq.1) then
            iiamin=1
         else
            iiamin=0
         end if
         return
      end if

      itmp=1
      imin=abs(iarray(1))
      if (incx.eq.1) then
         do i=2,n
            if (abs(iarray(i)).lt.imin) then
               if (iarray(i).ne.0) then
                  itmp=i
                  imin=abs(iarray(i))
               else
                  iiamin=i
                  return
               end if
            end if
         end do
      else
         ix=1
         do i=2,n
            ix=ix+incx
            if (abs(iarray(ix)).lt.imin) then
               if (iarray(ix).ne.0) then
                  itmp=i
                  imin=abs(iarray(ix))
               else
                  iiamin=i
                  return
               end if
            end if
         end do
      end if
      iiamin=itmp

      return
      end

