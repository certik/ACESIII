
c dx(i) *= da

      subroutine xscal(n,da,dx,incx)
      implicit none
      double precision da,dx(*)
      integer i,incx,p,n,nincx

      if (n.le.0 .or. incx.le.0 .or. da.eq.1.d0) return

      if (incx.eq.1) then

#ifdef _UNICOS
         do i = 1, n
            dx(i) = da*dx(i)
         end do
#else
         p = iand(n,3)
         do i = 1, p
            dx(i) = da*dx(i)
         end do
         do i = 1+p, n, 4
            dx(i)   = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
         end do
#endif

      else

         nincx = n*incx
         do i = 1, nincx, incx
            dx(i) = da*dx(i)
         end do

      end if

      return
c     end subroutine xscal
      end

