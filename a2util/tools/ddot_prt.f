c
c MIMICS:
c    double = factor*ddot()
c OUTPUT:
c    label dx(i) * dy(i) * factor => ddot
c
      subroutine ddot_prt(n,dx,incx,dy,incy,
     &                    label,factor)
      double precision dx(*),dy(*),dtmp
      integer i,incx,incy,n
      character*(*) label
      double precision factor
      if (n.le.0) return
      dtmp=0.0D0
      do i = 0,(n-1)
         dtmp=dtmp+dx(1+(i*incx))*dy(1+(i*incy))
         print '(A,4(A,F10.6))',
     &   label,' ',dx(1+(i*incx)),' * ',dy(1+(i*incy)),' * ',factor,
     &   ' => ',dtmp
      end do
      return
      end
