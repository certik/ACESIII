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
       subroutine cwork6642_unopt(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                         ne1,ne2,nf1,nf2,
     *                      a1, a2, b1, b2, c1, c2, d1, d2,
     *                      e1, e2, f1, f2,
     *                      inda, indb, indc, indd, inde, indf,
c
     *                      x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *                          nk1,nk2,nl1,nl2,
     *                          m1, m2, n1, n2, r1, r2, s1, s2,
     *                          k1, k2, l1, l2,
     *                      indx1,
c
     *                      x2,ni1,ni2,nj1,nj2,nu1,nu2,nv1,nv2,
     *                          i1, i2, j1, j2, u1, u2, v1, v2,
     *                      indx2,
     *                      cind, flopcount, scr1, scr2, scr3)
c-------------------------------------------------------------------------
c   Performs a "6642" contraction:
c      6 index output array  : y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,ne1:ne2,nf1:nf2)
c      6 index operand array : x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,nk1:nk2,nl1:nl2)
c      4 index operand array : x2(ni1:ni2,nj1:nj2,nr1:nr2,ns1:ns2)
c      2 index contraction.
c
c   I. e., contract two of the the 4 indices of the 2nd operand array out
c   of the first operand array, replacing the indices by the 2nd operand's
c   non-contracted indices.

c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *        ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,m1,m2,n1,n2,
     *        r1,r2,s1,s2,i1,i2,j1,j2,k1,k2,l1,l2
      integer nu1,nu2,nv1,nv2, u1, u2, v1, v2
      integer inda, indb, indc, indd, inde, indf, indx1(6), indx2(4), 
     *        indi, indj
      integer cind(2), flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,
     *                   ne1:ne2,nf1:nf2)
      double precision x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,
     *                    nk1:nk2,nl1:nl2)
      double precision x2(ni1:ni2,nj1:nj2,nu1:nu2,nv1:nv2)
      double precision scr1(*), scr2(*), scr3(*)

      integer ia,ib,ic,id,ie,if, ix_con1, ix_con2, ix(6), ix2(4), iy
      integer ia2, ib2, ic2, id2, ie2, if2 
      integer ix2_con1, ix2_con2
      integer bcon1, bcon2, econ1, econ2
      integer i, j, a, b, c, d, e, f, g, h
      integer lda, ldb, ldc, m, n, k, l
      integer next
      double precision x2val, xval, etemp 
      
      flopcount = 0

      do a = a1, a2
      do b = b1, b2
      do c = c1, c2
      do d = d1, d2
      do e = e1, e2
      do f = f1, f2
         y(a,b,c,d,e,f) = 0.
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c     return 

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y and x2 
c   indices.
c---------------------------------------------------------------------------

      ia = 0
      ib = 0
      ic = 0
      id = 0
      ie = 0
      if = 0

      do i = 1, 6
         if (indx1(i) .eq. inda) then 
            ia = i 
         else if (indx1(i) .eq. indb) then
            ib = i 
         else if (indx1(i) .eq. indc) then
            ic = i
         else if (indx1(i) .eq. indd) then
            id = i
         else if (indx1(i) .eq. inde) then
            ie = i
         else if (indx1(i) .eq. indf) then
            if = i
         else if (indx1(i) .eq. cind(1)) then
            ix_con1 = i
         else if (indx1(i) .eq. cind(2)) then
            ix_con2 = i
         else
            print *,'Error: Invalid index for x1 in cwork6642'
            print *,'X1 index is ',indx1(i),' y indices: ',inda, indb,
     *              indc, indd, inde, indf,
     *              ' X1 indices: ',(indx1(j),j=1,6),
     *              ' cind = ',cind(1),cind(2) 
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
      ie2 = 0
      if2 = 0

      do i = 1, 4
         if (indx2(i) .eq. inda) then 
            ia2 = i 
         else if (indx2(i) .eq. indb) then
            ib2 = i 
         else if (indx2(i) .eq. indc) then
            ic2 = i
         else if (indx2(i) .eq. indd) then
            id2 = i
         else if (indx2(i) .eq. inde) then
            ie2 = i
         else if (indx2(i) .eq. indf) then
            if2 = i
         else if (indx2(i) .eq. cind(1)) then
            ix2_con1 = i
         else if (indx2(i) .eq. cind(2)) then
            ix2_con2 = i
         else
            print *,'Error: Invalid index for x2 in cwork6642'
            print *,'X2 index is ',indx2(i),' y indices: ',inda, indb,
     *              indc, indd, inde, indf,
     *              ' X2 indices: ',(indx2(j),j=1,4),
     *              ' cind = ',cind(1), cind(2)
            call abort_job() 
         endif
      enddo

      if (cind(1) .eq. indx1(1)) then
         bcon1 = m1
         econ1 = m2
      else if (cind(1) .eq. indx1(2)) then
         bcon1 = n1
         econ1 = n2
      else if (cind(1) .eq. indx1(3)) then
         bcon1 = r1
         econ1 = r2
      else if (cind(1) .eq. indx1(4)) then
         bcon1 = s1
         econ1 = s2
      else if (cind(1) .eq. indx1(5)) then
         bcon1 = k1
         econ1 = k2
      else if (cind(1) .eq. indx1(6)) then
         bcon1 = l1
         econ1 = l2
      endif

      if (cind(2) .eq. indx1(1)) then
         bcon2 = m1
         econ2 = m2
      else if (cind(2) .eq. indx1(2)) then
         bcon2 = n1
         econ2 = n2
      else if (cind(2) .eq. indx1(3)) then
         bcon2 = r1
         econ2 = r2
      else if (cind(2) .eq. indx1(4)) then
         bcon2 = s1
         econ2 = s2
      else if (cind(2) .eq. indx1(5)) then
         bcon2 = k1
         econ2 = k2
      else if (cind(2) .eq. indx1(6)) then
         bcon2 = l1
         econ2 = l2
      endif

c     print *,'cwork6642: Unoptimized version lineno ',current_line
      do a = a1, a2
         if (ia  .ne. 0) ix(ia)   = a
         if (ia2 .ne. 0) ix2(ia2) = a
      do b = b1, b2
         if (ib  .ne. 0) ix(ib)   = b
         if (ib2 .ne. 0) ix2(ib2) = b
      do c = c1, c2
         if (ic  .ne. 0) ix(ic)   = c
         if (ic2 .ne. 0) ix2(ic2) = c
      do d = d1, d2
         if (id  .ne. 0) ix(id)   = d
         if (id2 .ne. 0) ix2(id2) = d
      do e = e1, e2
         if (ie  .ne. 0) ix(ie)   = e
         if (ie2 .ne. 0) ix2(ie2) = e
      do f = f1, f2
         if (if  .ne. 0) ix(if)   = f
         if (if2 .ne. 0) ix2(if2) = f

         y(a,b,c,d,e,f) = 0.
         etemp = 0.0 

         do i = bcon1, econ1 
            ix(ix_con1) = i
            ix2(ix2_con1) = i
         do j = bcon2, econ2 
            ix(ix_con2) = j
            ix2(ix2_con2) = j
c           y(a,b,c,d,e,f) = y(a,b,c,d,e,f) + 
            etemp = etemp + 
     *             x2(ix2(1),ix2(2),ix2(3),ix2(4))*
     *             x1(ix(1),ix(2),ix(3),ix(4),ix(5),ix(6))
         enddo
         enddo

         y(a,b,c,d,e,f) = etemp  
         
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
