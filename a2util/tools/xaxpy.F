
c dy(i) += da*dx(i)

      subroutine xaxpy(n,da,dx,incx,dy,incy)
      implicit none
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,p,n

      if (n.le.0 .or. da.eq.0.d0) return

      if (incx.eq.1.and.incy.eq.1) then

#ifdef _UNICOS
         do i = 1, n
            dy(i) = dy(i) + da*dx(i)
         end do
#else
         p = iand(n,3)
         do i = 1, p
            dy(i) = dy(i) + da*dx(i)
         end do
         do i = 1+p, n, 4
            dy(i)   = dy(i)   + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         end do
#endif

      else

         ix = 1
         iy = 1
         if (incx.lt.0) ix = 1 + (-n+1)*incx
         if (incy.lt.0) iy = 1 + (-n+1)*incy
         do i = 1, n
            dy(iy) = dy(iy) + da*dx(ix)
            ix = ix + incx
            iy = iy + incy
         end do

      end if

      return
c     end subroutine xaxpy
      end

