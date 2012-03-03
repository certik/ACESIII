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
      subroutine cwork2221(y,na1,na2,nb1,nb2,a1,a2,b1,b2,inda, indb, 
     *                     x1,ne1,ne2,nf1,nf2,e1,e2,f1,f2, 
     *                     indx1, x2,ni1,ni2,nj1,nj2,
     *                     i1,i2,j1,j2, indi, indj, scr1, scr2)
c-------------------------------------------------------------------------
c   Performs a "2221" contraction: 
c      2 index output array
c      2 index operand array
c      2 index operand array
c      1 index contraction.
c
c   I. e., contract one of the the 2 indices of the 2nd operand array out 
c   of the first operand array, replacing the index by the 2nd operand's 
c   2nd index.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      integer na1,na2,nb1,nb2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,e1,e2,f1,f2,
     *        i1,i2,j1,j2
      integer inda, indb, indx1(2), indi, indj

      double precision y(na1:na2,nb1:nb2)
      double precision x1(ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2)
      double precision scr1(*), scr2(*)

      integer i, j, a, b, e, f
#ifdef __crayx1
      integer*8 m, n, k, lda, ldb, ldc
#else
      integer m, n, k, lda, ldb, ldc
#endif
      integer next
      character*1 op_a, op_b

      next = 0
      do f = f1,f2
      do e = e1,e2
         next = next + 1
         scr1(next) = x1(e,f)
      enddo
      enddo

      next = 0
      do j = j1,j2
      do i = i1,i2
         next = next + 1
         scr2(next) = x2(i,j)
      enddo
      enddo

      m = a2-a1+1
      n = b2-b1+1
      lda = e2-e1+1
      ldc = na2-na1+1
      ldb = i2-i1+1

      op_a = 'X'
      op_b = 'X'
      if (inda .eq. indx1(1)) then
         k = f2-f1+1
         if (indb .eq. indi) then

c-------------------------------------------------------------------------
c   Y(a,b) = x1(a,c) * x2(b,c)
c-------------------------------------------------------------------------

            op_a = 'N'
            op_b = 'T'
         else

c-------------------------------------------------------------------------
c   Y(a,b) = x1(a,c) * x2(c,b)
c-------------------------------------------------------------------------

            op_a = 'N'
            op_b = 'N'
         endif
      else
         k = e2-e1+1
         if (indb .eq. indi) then

c-------------------------------------------------------------------------
c   Y(a,b) = x1(c,a)*x2(b,c)
c-------------------------------------------------------------------------

            op_a = 'T'
            op_b = 'T'
         else
                                                                                
c-------------------------------------------------------------------------
c   Y(a,b) = x1(c,a)*x2(c,b)
c-------------------------------------------------------------------------

            op_a = 'T'
            op_b = 'N'
         endif
      endif

      if (op_a .eq. 'X' .or. op_b .eq. 'X') then
         print *,'Error in cwork2221 logic, line ',current_line
         call abort_job()
      endif

      call dgemm(op_a, op_b, m,n,k,1.d0, scr1, lda,
     *                 scr2, ldb,
     *                 0.d0, y(a1,b1), ldc)
      return
      end
