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
      subroutine flowqsdmin(scratch, grad, ehess, eigvh, dmr1, dmr2,
     &                      stp, eps, nx, ndim, ncycle)
c
c This subroutine implements the quadratic gradient descent
c algorithm of J.-Q. Sun and K.Ruedenberg in JCP 99, 5257 (1993).
c
      implicit double precision(a-h, o-z)
#include "fnamelen.par"
      character*(fnamelen) fname
      logical yesno, update, nwtonrpson, negeval, genstpsz
c
      dimension scratch(nx*nx), grad(ndim), ehess(ndim, ndim), 
     &          eigvh(ndim, ndim), dmr1(ndim, ndim), dmr2(ndim,ndim)
c
c find the eigenvalues and eigenvectors of the hessian to build the
c solution to eqn. 11 in Sun & Ruedenberg.
c
      call eig(ehess, eigvh, ndim, ndim, 1)
c
c Find the tangential vector, the scalar curvature, and define some
c ad hoc factors from Sun & Ruedenburg.
c
      update     = .false.
      nwtonrpson = .false.
      negeval    = .false.
      genstpsz    = .true.
      unew    = 0.0d0
      u_upper = 0.0d0
      u_lower = 0.0d0
      hvalue  = 0.0d0
      anu     = 5.0d0/6.0d0
c
      call qsdline(scratch, ehess, eigvh, dmr1, dmr2, unew, u_upper,
     &             u_lower, hvalue, curvtre, ndim, nx, update,
     &             nwtonrpson, negeval, genstpsz)
c
c construct the curvature dependent step size (eqn. 20b)
c
      crvdstsz = 0.25d0 + (0.75d0)/(1.d0 + 3.0d0*curvtre**anu)
      Write(6, *) "The value of curv =", curvtre
      eh = 0.3d0*crvdstsz
      Write(6, *) "The value of stpsize = ", eh
c      
c
c       if (ncycle .eq. 1) then
c
c P_0=Q_0=R_0 at the initial point. Calculate P* and Q_1 and 
c Save P*
c
c           call xdcopy(ndim, scratch(1), 1, scratch(1 + 2*ndim), 1)
           
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                  (eh/2)**2, eps, nx, ndim)
c
c          open(unit = 99, file = "QSDPINTS", form = "formatted")
c          write(99, 999) eh
c          write(99, 9999) (scratch(i + 2*ndim), i = 1, ndim)
c
c Get the Q_1 based on P*
c
          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
     &                   (eh/2.0d0), eps, nx, ndim)

c          call vadd(scratch(1 + 2*ndim), scratch(1 + 2*ndim), 
c     &              scratch(1), ndim, -1.0d0)
c
c          call xdcopy(ndim, scratch(1 + 2*ndim), 1, scratch(1 + ndim),
c     &               1)
c
c       else if (ncycle .eq. 2) then
c
c Calculate R_1 based on P*
c             
c          call gfname('QSDPINTS', fname, ilength)
c          inquire(file=fname(1:ilength), exist=yesno)
c          if (yesno) open(unit = 99, file = "QSDPINTS", 
c     &                    form = "formatted", status ="old")
c
c          read(99, 999) ehp
c          read(99, 9999) (scratch(i + 2*ndim), i = 1, ndim)
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                   (ehp/4.0d0)**2, eps, nx, ndim)
c
c Calculate P_1 based on Q_1
c
c          call xdcopy(ndim, scratch(1), 1, scratch(1 + 2*ndim), 1)
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                  (ehp/2.0d0)**2, eps, nx, ndim)
c          rewind(99)
c          write(99, 999) eh
c          write(99, 9999) (scratch(i + 2*ndim), i = 1, ndim)
c
c Calculate Q_2 based on P_1
c
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                  (eh/2.0d0)**2, eps, nx, ndim)
c
c          call vadd(scratch(1 + 2*ndim), scratch(1 + 2*ndim), 
c     &             scratch(1), ndim, -1.0d0)
c
c          call xdcopy(ndim, scratch(1 + 2*ndim), 1, scratch(1 + ndim),
c     &               1)
c
c       else
c
c Calculate R_n based on P_n-1
c
c          call gfname('QSDPINTS', fname, ilength)
c          inquire(file=fname(1:ilength), exist=yesno)
c          if (yesno) open(unit = 99, file = "QSDPINTS", 
c     &                    form = "formatted", status="old")
c
c          read(99, 999) ehp
c          read(99, 9999) (scratch(i + 2*ndim), i = 1, ndim)
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                  (ehp/2.0d0)**2, eps, nx, ndim)
c
c Calculate P_n based on Q_n
c
c          call xdcopy(ndim, scratch(1), 1, scratch(1 + 2*ndim), 1)
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                  (eh/2.0d0)**2, eps, nx, ndim)
c          rewind(99)
c          write(99, 999) eh
c          write(99, 9999) (scratch(i + 2*ndim), i = 1, ndim)
c
c Calculate Q_n+1 based on P_n
c
c          call getptspqr(scratch, ehess, eigvh, dmr1, dmr2, 
c     &                  (eh/2.0d0)**2, eps, nx, ndim)
c
c           call vadd(scratch(1 + 2*ndim), scratch(1 + 2*ndim), 
c     &             scratch(1), ndim, -1.0d0)
c
c           call xdcopy(ndim, scratch(1 + 2*ndim), 1, 
c     &                scratch(1 + ndim), 1)
c
c
c      endif
c       
 999    format(1x, f15.8)
 9999   format(1x, 3(f15.8, 1x))
c        
       return
       end
