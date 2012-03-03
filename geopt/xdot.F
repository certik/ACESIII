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

