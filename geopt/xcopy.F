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

#include "f_types.h"
c#define _ABORT_ON_OVERLAP

      subroutine xcopy(n,dx,incx,dy,incy)
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,n
      F_ADR c_adr, z
      external c_adr
      if(n.le.0)return
      if ((incx.eq.1).and.(incy.eq.1)) then
         if (c_adr(dy).eq.c_adr(dx)) return
#ifdef _DEBUG
c      o complain about overlapping data
#ifdef _PTRS_ARE_WORDS
c        This might be a defect. I am assuming c_adr returns 64-bit words.
c        I could divide by iintfp, but that would require including /MACHSP/,
c        which means xcopy would have to be moved out of the tools library.
         z = abs(c_adr(dx)-c_adr(dy))
#else
         z = abs(c_adr(dx)-c_adr(dy))/8
#endif
         if (z.lt.n) then
            print '(/)'
            print *, '@XCOPY: WARNING - the source and destination ',
     &               'data overlap'
            print *, '        src address: ',c_adr(dx)
            print *, '        dst address: ',c_adr(dy)
            print *, '        difference:  ',z,' doubles'
            print *, '        data length: ',n,' doubles'
            print '(/)'
#ifdef _ABORT_ON_OVERLAP
            call c_abort
#endif
         end if
#endif
         call c_memmove(dy,dx,n*8)
      else
         ix = 1
         iy = 1
         if(incx.lt.0)ix = (-n+1)*incx + 1
         if(incy.lt.0)iy = (-n+1)*incy + 1
         do i = 1,n
           dy(iy) = dx(ix)
           ix = ix + incx
           iy = iy + incy
         end do
      end if
      return
      end

