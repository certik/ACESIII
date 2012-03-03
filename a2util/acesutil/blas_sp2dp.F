
#ifndef USE_SP_BLAS
      subroutine saxpy(n,da,dx,incx,dy,incy)
      double precision dx(*),dy(*),da
      integer incx,incy,n
      call daxpy(n,da,dx,incx,dy,incy)
      return
      end

      subroutine scopy(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
      integer incx,incy,n
      call dcopy(n,dx,incx,dy,incy)
      return
      end

      subroutine sscal(n,da,dx,incx)
      double precision da,dx(*)
      integer incx,n
      call dscal(n,da,dx,incx)
      return
      end

      double precision function sdot(n,dx,incx,dy,incy)
      double precision dx(*),dy(*),ddot
      integer incx,incy,n
      sdot=ddot(n,dx,incx,dy,incy)
      return
      end

      double precision function snrm2(n,x,incx)
      double precision x(*),dnrm2
      integer incx,n
      snrm2=dnrm2(n,x,incx)
      return
      end

      double precision function ssum(n,dx,incx)
      double precision dx(*),dsum
      integer incx,n
      ssum=dsum(n,dx,incx)
      return
      end

      integer function isamax(n,dx,incx)
      double precision dx(*)
      integer incx,n,idamax
      isamax=idamax(n,dx,incx)
      return
      end
#endif /* !USE_SP_BLAS */
