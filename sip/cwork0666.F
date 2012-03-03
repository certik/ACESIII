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
      subroutine cwork0666(y, x1,
     *                     ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *                     nm1,nm2,nn1,nn2, 
     *                     e1,e2,f1,f2,g1,g2,h1,h2,
     *                     m1,m2,n1,n2, 
     *                     indx1, x2,
     *                     ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                     no1,no2,np1,np2, 
     *                     i1,i2,j1,j2,k1,k2,l1,l2,i5,j5,i6,j6,
     *                     indx2, flopcount)
c-------------------------------------------------------------------------
c   Performs a "0666" contraction: 
c      0 index output array
c      6 index operand array
c      6 index operand array
c      6 index contraction.
c
c   I. e., contract all 4 indices of the 2nd operand array with the
c   first operand array, storing the result in a scalar. 
c--------------------------------------------------------------------------
      implicit none

      integer ne1,ne2,nf1,nf2,
     *        ng1,ng2,nh1,nh2,
     *        ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2, 
     *        nm1,nm2,nn1,nn2,
     *        no1,no2,np1,np2 
      integer e1,e2,f1,f2,g1,g2,h1,h2,
     *        i1,i2,j1,j2,k1,k2,l1,l2, 
     *        m1,m2,n1,n2,i5,j5,i6,j6  
      integer indx1(4), indx2(4), indi, indj
      integer flopcount

      double precision xddot
      double precision y
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2,
     *                    nm1:nm2,nn1:nn2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2,
     *                    no1:no2,np1:np2)

      integer ii,ij,ik,il
      integer i, j, k, l, m, n, ix(4)
      integer one 

      flopcount = 0
      y = 0.0d0 

c--------------------------------------------------------------------------
c   Indices of x1, x2 line up exactly.  Use optimized contraction loop.
c--------------------------------------------------------------------------

         n = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(l2-l1+1)
     *      *(j5-i5+1)*(j6-i6+1)
         one = 1
         y = xddot(n, x1(i1,j1,k1,l1,i5,i6), one, 
     *                x2(i1,j1,k1,l1,i5,i6), one)
c     y = 0.0d0 
c     do i = i1, i2 
c     do j = j1, j2 
c     do k = k1, k2 
c     do l = l1, l2 
c     do m = i5, j5  
c     do n = i6, j6  
c        y = y + x1(i,j,k,l,m,n)*x2(i,j,k,l,m,n) 
c     enddo 
c     enddo 
c     enddo 
c     enddo 
c     enddo 
c     enddo 

      return
      end
