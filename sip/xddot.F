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
      double precision function xddot(n,dx,incx,dy,incy)

      double precision dx(*),dy(*),dsum
      integer incx, incy, n

      if (n.le.0) return

      dsum = 0.d0
      do i = 1, n
         dsum = dsum + dx(1+(i-1)*incx) * dy(1+(i-1)*incy)
      enddo

      xddot = dsum
      return
      end

