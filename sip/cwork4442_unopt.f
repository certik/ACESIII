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
       subroutine cwork4442_unopt(y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                     a1,a2,b1,b2,c1,c2,d1,d2,inda, indb, 
     *                     indc, indd, x1,
     *                     ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *                     e1,e2,f1,f2,g1,g2,h1,h2, 
     *                     indx1, x2,
     *                     ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                     i1,i2,j1,j2,k1,k2,l1,l2,
     *                     indx2, cind, flopcount, scr1, scr2) 
c-------------------------------------------------------------------------
c   Performs a "4442" contraction: 
c      4 index output array
c      4 index operand array
c      4 index operand array
c      2 index contraction.
c
c   I. e., contract two of the the 4 indices of the 2nd operand array out 
c   of the first operand array, replacing the indices by the 2nd operand's 
c   non-contracted indices.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ng1,ng2,nh1,nh2,
     *        ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,
     *        i1,i2,j1,j2,k1,k2,l1,l2
      integer inda, indb, indc, indd, indx1(4), indx2(4), indi, indj
      integer cind(2), flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
      double precision scr1(*), scr2(*)

      integer ia, ib, ic, id, ix_con1, ix_con2, ix(4), ix2(4), iy
      integer ia2, ib2, ic2, id2
      integer ix2_con1, ix2_con2
      integer bcon1, bcon2, econ1, econ2
      integer i, j, a, b, c, d, e, f, g, h
      integer lda, ldb, ldc, m, n, k, l
      integer next
      double precision x2val, xval
      
      flopcount = 0

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y and x2 
c   indices.
c---------------------------------------------------------------------------

      ia = 0
      ib = 0
      ic = 0
      id = 0

      do i = 1, 4
         if (indx1(i) .eq. inda) then 
            ia = i 
         else if (indx1(i) .eq. indb) then
            ib = i 
         else if (indx1(i) .eq. indc) then
            ic = i
         else if (indx1(i) .eq. indd) then
            id = i
         else if (indx1(i) .eq. cind(1)) then
            ix_con1 = i
         else if (indx1(i) .eq. cind(2)) then
            ix_con2 = i 
         else
            print *,'Error: Invalid index for x1 in cwork4442'
            print *,'X1 index is ',indx1(i),' y indices: ',inda, indb,
     *              indc, indd,
     *              ' X1 indices: ',(indx1(j),j=1,4),
     *              ' cind = ',cind(1), cind(2) 
            call abort_job() 
         endif
      enddo

c---------------------------------------------------------------------------
c   Find which indices of the "x2" operand match the various y and x1
c   indices.
c---------------------------------------------------------------------------

      ia2 = 0
      ib2 = 0
      ic2 = 0
      id2 = 0

      do i = 1, 4
         if (indx2(i) .eq. inda) then 
            ia2 = i 
         else if (indx2(i) .eq. indb) then
            ib2 = i 
         else if (indx2(i) .eq. indc) then
            ic2 = i
         else if (indx2(i) .eq. indd) then
            id2 = i
         else if (indx2(i) .eq. cind(1)) then
            ix2_con1 = i
         else if (indx2(i) .eq. cind(2)) then
            ix2_con2 = i
         else
            print *,'Error: Invalid index for x2 in cwork4442'
            print *,'X2 index is ',indx2(i),' y indices: ',inda, indb,
     *              indc, indd,
     *              ' X2 indices: ',(indx2(j),j=1,4),
     *              ' cind = ',cind(1),cind(2)
            call abort_job() 
         endif
      enddo

      if (cind(1) .eq. indx1(1)) then
         bcon1 = e1
         econ1 = e2
      else if (cind(1) .eq. indx1(2)) then
         bcon1 = f1
         econ1 = f2
      else if (cind(1) .eq. indx1(3)) then
         bcon1 = g1
         econ1 = g2
      else if (cind(1) .eq. indx1(4)) then
         bcon1 = h1
         econ1 = h2
      endif

      if (cind(2) .eq. indx1(1)) then
         bcon2 = e1
         econ2 = e2
      else if (cind(2) .eq. indx1(2)) then
         bcon2 = f1
         econ2 = f2
      else if (cind(2) .eq. indx1(3)) then
         bcon2 = g1
         econ2 = g2
      else if (cind(2) .eq. indx1(4)) then
         bcon2 = h1
         econ2 = h2
      endif

      print *,'cwork4442: Unoptimized version lineno ',current_line
      do a = a1, a2
         if (ia .ne. 0) ix(ia)  = a
         if (ia2 .ne. 0) ix2(ia2) = a
      do b = b1, b2
         if (ib .ne. 0) ix(ib) = b
         if (ib2 .ne. 0) ix2(ib2) = b
      do c = c1, c2
         if (ic .ne. 0) ix(ic) = c
         if (ic2 .ne. 0) ix2(ic2) = c
      do d = d1, d2
         if (id .ne. 0) ix(id) = d
         if (id2 .ne. 0) ix2(id2) = d
         y(a,b,c,d) = 0.

         do j = bcon2, econ2
            ix(ix_con2) = j 
            ix2(ix2_con2) = j 
         do i = bcon1, econ1 
            ix(ix_con1) = i
            ix2(ix2_con1) = i
            y(a,b,c,d) = y(a,b,c,d) + 
     *             x2(ix2(1),ix2(2),ix2(3),ix2(4))*
     *             x1(ix(1),ix(2),ix(3),ix(4))
         enddo
         enddo
         
      enddo
      enddo
      enddo
      enddo

      return
      end
