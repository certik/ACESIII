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
      subroutine qsdline(scratch, ehess, eigvh, dmr1, dmr2, unew,
     &                   u_upper, u_lower, hvalue, curvtre, ndim, nx,
     &                   update, nwtonrpson, negeval, genstpsz)
c
c this routine calculate quadaratic steepest descent line (eqn. 11
c sun & rudernberg jcp, 99, 5258, 1993.
c
      implicit double precision (a-h, o-z)
      logical update, nwtonrpson, negeval, genstpsz
c
      dimension scratch(nx*nx), eigvh(ndim, ndim), dmr1(ndim, ndim),
     &          dmr2(ndim, ndim), ehess(ndim, ndim)
c
      call xgemm('t','n', ndim, 1, ndim, 1.0d0, eigvh, ndim,
     &            scratch(1 + ndim), ndim, 0.0d0,
     &            scratch(1 + 2*ndim), ndim)
c
      if (genstpsz) then
         call curve(scratch, ehess, curvtre, ndim, nx)
         return
      endif

      if (negeval) then
         u_lower =  dexp(-dlog(1.0d0 + 2.0d0*hvalue/
     &              dabs(scratch(1 + 2*ndim)))/dabs(ehess(1, 1)))
         unew = u_lower
      endif
c
      do 10 i = 1, ndim
         dmr1(i, i) = unew**ehess(i, i)
 10   continue
c
      do 80 i = 1, ndim
         scratch(i + 3*ndim) = dmr1(i, i)*scratch(i + 2*ndim)
  80   continue
c
       if (negeval) then
          call xdcopy(ndim,scratch(1 + 3*ndim),1,scratch(1 + 4*ndim),
     &               1)
          return
       endif
c
c Let's first check the search can be reduced to the Newton-Rapshon search
c
      dtmp = xdot(ndim,scratch(1 + 2*ndim),1,scratch(1 + 2*ndim),1)
      threshld = dsqrt(dtmp)

      If (threshld .lt. hvalue) then
         write(6, *) "  QSD search reduce to Newton-Raphson update"
         call vminus(scratch(1+ ndim), ndim)
         nwtonrpson = .true.
         return
      endif

      if (.not. update) then
c
         call vadd(scratch(1 + 3*ndim), scratch(1 + 3*ndim),
     &             scratch(1 + 2*ndim), ndim, -1.0d0)

         dtmp = xdot(ndim,scratch(1 + 3*ndim),1,scratch(1 + 3*ndim),1)
         xupqr = dsqrt(dtmp)
c
         if (xupqr .gt. hvalue) then
            u_lower = unew
            call xdcopy(ndim,scratch(1 + 3*ndim),1,scratch(1 + 4*ndim),
     &                 1)
         else
            u_upper = unew
            call xdcopy(ndim,scratch(1 + 3*ndim),1,scratch(1 + 5*ndim),
     &                 1)
         endif
c
      else
c
         call xgemm('n','n', ndim, 1, ndim, 1.0d0, eigvh, ndim,
     &               scratch(1 + 3*ndim), ndim, 0.d0,
     &               scratch(1 + 2*ndim), ndim)
c
         call vadd(scratch(1 + 3*ndim), scratch(1 + 2*ndim),
     &             scratch(1 + ndim), ndim, -1.0d0)

         call xdcopy(ndim, scratch(1 + 3*ndim), 1, scratch(1 + ndim), 1)

      endif
c
      return
      end
