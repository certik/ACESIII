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

c complex division, (cr,ci) = (ar,ai)/(br,bi)

      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
      double precision s,ars,ais,brs,bis,s_inv
c old
c      s = dabs(br) + dabs(bi)
c      ars = ar/s
c      ais = ai/s
c      brs = br/s
c      bis = bi/s
c      s = brs**2 + bis**2
c      cr = (ars*brs + ais*bis)/s
c      ci = (ais*brs - ars*bis)/s
c new
      s_inv = 1.0d0 / ( dabs(br) + dabs(bi) )
      ars = ar * s_inv
      ais = ai * s_inv
      brs = br * s_inv
      bis = bi * s_inv
      s_inv = 1.0d0 / ( brs*brs + bis*bis )
      cr = (ars*brs + ais*bis) * s_inv
      ci = (ais*brs - ars*bis) * s_inv
c end
      return
      end
