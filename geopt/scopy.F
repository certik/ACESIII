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
      subroutine scopy(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
#ifdef ACES64
      integer*8 n, incx, incy
#else
      integer n, incx, incy
#endif
      integer i,ix,iy

         ix = 1
         iy = 1
         if(incx.lt.0)ix = (-n+1)*incx + 1
         if(incy.lt.0)iy = (-n+1)*incy + 1
         do i = 1,n
           dy(iy) = dx(ix)
           ix = ix + incx
           iy = iy + incy
         end do
      return
      end

