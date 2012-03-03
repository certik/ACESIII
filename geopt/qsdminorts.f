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
      subroutine qsdminorts(scratch, grad, ehess, eigvh, dmr1, dmr2,
     &                      stp, eps, igts, qsd, nx, ndim, ncycle)
c
c This subroutine implements the quadratic image gradient and hessian
c algorithm of J.-Q. Sun and K.Ruedenberg in JCP 101, 2157 (1994).
c
c Input variables: 
c ndim - Integer number of coordinates
c grad - Real vector of gradients
c ehess - Real matrix of hessian elements
c dmr# - Dummy real matrix
c scratch - Dummy real vector
c
c  Output variables:
c  grad - Real vector of image gradient
c  ehess - Real matrix of image hessian
c
c  Date written : 12/14/99
c
c  Last modified :
c
c  Author : Keith Runge, Ajith Perera
c  Modifications and adaptation to ACES II 04/03/2000
c
      implicit double precision (a-h, o-z)
c
      logical qsd, igts
c
      dimension scratch(nx*nx), grad(ndim), ehess(ndim, ndim),
     &           eigvh(ndim, ndim), dmr1(ndim, ndim), dmr2(ndim, ndim)
c
c Rebuild the hessian in symmetry coordinates.
c
      do 10 i = 1, ndim
         call xgemm('n', 't', ndim, ndim, 1, ehess(i, i), eigvh(1, i),
     &               ndim, eigvh(1, i), ndim, 1.0d0, dmr1,  ndim)
 10   continue
c
      call xdcopy(ndim*ndim, dmr1, 1, ehess, 1)
c
c Find the bond critical point (H^-1g)
c
      call minv(dmr1, ndim, ndim, dmr2, Det, 1.0d-8, 0, 1)
      call xgemm('n', 'n', ndim, 1, ndim, 1.d0, dmr1, ndim, grad, ndim,
     &            0.d0, scratch(ndim + 1), ndim)
c
      if (igts) then
c
c Build the image hessian and image gradient.
c    
         call zero(dmr2, ndim*ndim)

         do 20 i = 1, ndim
            dmr2(i, i) = 1.0d0
c            do 30 j = 1, ndim
c               dmr2(i, j) = dmr2(i, j) - 2.d0*eigvh(i,1)*eigvh(j,1)
c 30         continue
 20      continue

         call xgemm('n', 't', ndim, ndim, 1, -2.0d0, eigvh(1, 1),
     &               ndim, eigvh(1, 1), ndim, 1.0d0, dmr2,  ndim)
c
c The image hessian
c    
         call xgemm('n', 'n', ndim, ndim, ndim, 1.d0, ehess, ndim, dmr2,
     &               ndim, 0.d0, dmr1, ndim)
c
c The image gradient
c    
         call xgemm('n', 'n', ndim, 1, ndim, 1.d0, dmr1, ndim, grad,
     &               ndim, 0.d0, scratch(1 + 2*ndim), ndim)
c
         call xdcopy(ndim, scratch(1+ 2*ndim) , 1, grad, 1)
         call xdcopy(ndim*ndim, dmr1, 1, ehess, 1)
c
c Use the Quadratic Steepest Descent method to locate the image minimum.
c
         call flowqsdmin(scratch, grad, ehess, eigvh, dmr1, dmr2, stp,
     &                   eps, nx, ndim, ncycle)
c    
      else
c
c Look for a minimum using Quadratic Steepest Descent.
c
         call flowqsdmin(scratch, grad, ehess, eigvh, dmr1, dmr2, stp,
     &                   eps, nx, ndim, ncycle)
c        
      endif
c     
      return
      end
