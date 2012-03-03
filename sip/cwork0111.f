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
      subroutine cwork0111(y,
     *                     x1,ne1,ne2,e1,e2, indx1,
     *                     x2,ni1,ni2,i1,i2, indi, 
     *                     flopcount) 
c-------------------------------------------------------------------------
c   Performs a "0111" contraction: 
c      0 index output array (i. e. scalar)
c      1 index operand array
c      1 index operand array
c      1 index contraction.
c
c   I. e., contract the 1 index of the 2nd operand array out of the 
c   first operand array, forming a scalar result.
c--------------------------------------------------------------------------
      implicit none

      integer ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer e1,e2,f1,f2,i1,i2,j1,j2
      integer indx1(2), indi, indj
      integer flopcount

      double precision y
      double precision x1(ne1:ne2)
      double precision x2(ni1:ni2)

      integer ii, ij, ix(2)
      integer i, j, k, n

c--------------------------------------------------------------------------
c   Perform contraction. 
c--------------------------------------------------------------------------

      n = (i2-i1+1)*(j2-j1+1)

      y = 0.
      do i = i1, i2
         y = y + x1(i)*x2(i)
      enddo 

      return
      end
