c
c MIMICS:
c    dy(i) += (factor*da)*dx(i)
c OUTPUT:
c    label(start+(i*step)) += factor * da * dx(1+(i*incx))
c
      subroutine daxpy_prt(n,da,dx,incx,dy,incy,
     &                     label,start,step,factor)
      double precision dx(*),dy(*),da
      integer i,incx,incy,n
      character*(*) label
      integer start, step
      double precision factor
      if (n.le.0) return
      if (da.eq.(0.0d0)) return
      do i = 0,(n-1)
         print '(A,A,I3,3(A,F10.6))',
     &      label,'(',start+(i*step),') += ',
     &      factor,' * ',da,' * ',dx(1+(i*incx))
      end do
      return
      end
