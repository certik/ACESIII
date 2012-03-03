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
      subroutine cwork0333(y, x1,
     *                     ne1,ne2,nf1,nf2,ng1,ng2,
     *                     e1,e2,f1,f2,g1,g2,
     *                     indx1, x2,
     *                     ni1,ni2,nj1,nj2,nk1,nk2,
     *                     i1,i2,j1,j2,k1,k2,
     *                     indx2, flopcount)
c-------------------------------------------------------------------------
c   Performs a "0333" contraction: 
c      0 index output array
c      3 index operand array
c      3 index operand array
c      3 index contraction.
c
c   I. e., contract all 3 indices of the 2nd operand array with the
c   first operand array, storing the result in a scalar. 
c--------------------------------------------------------------------------
      implicit none

      integer ne1,ne2,nf1,nf2,ng1,ng2,
     *        ni1,ni2,nj1,nj2,nk1,nk2
      integer e1,e2,f1,f2,g1,g2,
     *        i1,i2,j1,j2,k1,k2
      integer indx1(4), indx2(4), indi, indj
      integer flopcount

      double precision xddot
      double precision y
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2)

      integer ii,ij,ik
      integer i, j, k, ix(3)
      integer one, n

      flopcount = 0

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various x2 
c   indices.
c---------------------------------------------------------------------------

      do i = 1, 3
         if (indx1(i) .eq. indx2(1)) ii = i
         if (indx1(i) .eq. indx2(2)) ij = i
         if (indx1(i) .eq. indx2(3)) ik = i
      enddo

      flopcount = 2*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)

      if (ii .ne. 1 .or.
     *    ij .ne. 2 .or.
     *    ik .ne. 3) then  

c--------------------------------------------------------------------------
c   Use mapping array to perform contraction via indirect addressing.
c--------------------------------------------------------------------------

         y = 0.

         do k = k1, k2
            ix(ik) = k
         do j = j1, j2
            ix(ij) = j 
         do i = i1, i2 
            ix(ii) = i
            y = y + x2(i,j,k)*x1(ix(1),ix(2),ix(3))
         enddo
         enddo
         enddo

      else

c--------------------------------------------------------------------------
c   Indices of x1, x2 line up exactly.  Use optimized contraction loop.
c--------------------------------------------------------------------------

         n = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
         one = 1
         y = xddot(n, x1(i1,j1,k1), one, x2(i1,j1,k1), one)
      endif

      return
      end
