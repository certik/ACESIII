
      double precision function xdot(n,dx,incx,dy,incy)

      double precision dx(*),dy(*),ddot
      integer incx, incy, n

#ifdef _32BIT_BLAS
      integer*4 int_n, int_incx, int_incy
#else
      integer   int_n, int_incx, int_incy
#endif /* _32BIT_BLAS */

      if (n.le.0) return

      int_incx = incx
      int_incy = incy
      int_n    = n

#ifdef USE_SP_BLAS
      xdot = sdot(int_n,dx,int_incx,dy,int_incy)
#else
      xdot = ddot(int_n,dx,int_incx,dy,int_incy)
#endif

      return
      end

