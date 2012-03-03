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
      subroutine daxpy(n,da,dx,incx,dy,incy)
      implicit none
      double precision dx(*),dy(*),da
#ifdef ACES64
      integer*8 i,incx,incy,n
#else
      integer i,incx,incy,n
#endif
      integer ix, iy

      if (n .le. 0 .or. da .eq. 0.d0) return

      iy = 1
      ix = 1
      if (incx.lt.0) ix = 1 + (-n+1)*incx
      if (incy.lt.0) iy = 1 + (-n+1)*incy
 
      do i = 1, n
         dy(iy) = dy(iy) + da * dx(ix)
         iy = iy + incy
         ix = ix + incx
      enddo

      return
      end
