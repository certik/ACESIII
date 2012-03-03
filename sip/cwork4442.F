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
       subroutine cwork4442(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                          a1, a2, b1, b2, c1, c2, d1, d2,
     *                      inda, indb, indc, indd,
     *                      x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *                          e1, e2, f1, f2, g1, g2, h1, h2, 
     *                      indx1,
     *                      x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                          i1, i2, j1, j2, k1, k2, l1, l2,
     *                      indx2,
     *                      cind, flopcount, scr1, scr2, scr3) 
c-------------------------------------------------------------------------
c   Performs a "4442" contraction: 
c      4 index output array  : y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
c      4 index operand array : x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2)
c      4 index operand array : x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
c      2 index contraction.
c
c   I. e., contract two of the the 4 indices of the 2nd operand array out 
c   of the first operand array, replacing the indices by the 2nd operand's 
c   non-contracted indices.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *        ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *        ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2
      integer a1,a2,b1,b2,c1,c2,d1,d2,
     *        e1,e2,f1,f2,g1,g2,h1,h2,
     *        i1,i2,j1,j2,k1,k2,l1,l2
      integer inda, indb, indc, indd, indx1(4), indx2(4)
      integer cind(2), flopcount

c      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision y(*)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
      double precision scr1(*), scr2(*), scr3(*)

      integer p1, p2, q1, q2
      integer na, nb, nc, nd, ne, nf, ng, nh, ni, nj, nk, nl, np, nq

      integer i1a, i1b, i1c, i1d, i1p, i1q
      integer i2a, i2b, i2c, i2d, i2p, i2q
      integer a, b, c, d, e, f, g, h, i, j, l, p, q

#ifdef __crayx1
      integer*8 lda, ldb, ldc, m, nn, k
#else
      integer lda, ldb, ldc, m, nn, k
#endif
      integer n
      character*1 op_a, op_b
      logical b2pq
      logical debug

      debug = cwork_debug
      flopcount = 0

      if (.false.) then

c   o assert x1 is contiguous
      if (e1.ne.ne1.or.f1.ne.nf1.or.g1.ne.ng1.or.h1.ne.nh1.or.
     &    e2.ne.ne2.or.f2.ne.nf2.or.g2.ne.ng2.or.h2.ne.nh2    ) then
         print *, '@CWORK4442: Assertion failed.'
         print *, '            The x1 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert x2 is contiguous
      if (i1.ne.ni1.or.j1.ne.nj1.or.k1.ne.nk1.or.l1.ne.nl1.or.
     &    i2.ne.ni2.or.j2.ne.nj2.or.k2.ne.nk2.or.l2.ne.nl2    ) then
         print *, '@CWORK4442: Assertion failed.'
         print *, '            The x2 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert y is contiguous
      if (a1.ne.na1.or.b1.ne.nb1.or.c1.ne.nc1.or.d1.ne.nd1.or.
     &    a2.ne.na2.or.b2.ne.nb2.or.c2.ne.nc2.or.d2.ne.nd2    ) then
         print *, '@CWORK4442: Assertion failed.'
         print *, '            The product tensor is not contiguous.'
         call abort_job()
      end if

c     end if (.false.)
      end if

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y and x2 
c   indices.
c---------------------------------------------------------------------------

      i1a = 0
      i1b = 0
      i1c = 0
      i1d = 0
      i1p = 0
      i1q = 0

      do i = 1, 4
         if      (indx1(i) .eq. inda) then 
            i1a = i 
         else if (indx1(i) .eq. indb) then
            i1b = i 
         else if (indx1(i) .eq. indc) then
            i1c = i
         else if (indx1(i) .eq. indd) then
            i1d = i
         else if (indx1(i) .eq. cind(1)) then
            i1p = i
         else if (indx1(i) .eq. cind(2)) then
            i1q = i 
         else
            print *,'Error: Invalid index for x1 in cwork4442'
            print *,'X1 index is ',indx1(i),
     *              ' y indices: ',inda,indb,indc,indd,
     *              ' X1 indices: ',(indx1(j),j=1,4),
     *              ' cind = ',cind(1), cind(2) 
            call abort_job() 
         end if
      end do

c---------------------------------------------------------------------------
c   Find which indices of the "x2" operand match the various y and x1
c   indices.
c---------------------------------------------------------------------------

      i2a = 0
      i2b = 0
      i2c = 0
      i2d = 0
      i2p = 0
      i2q = 0

      do i = 1, 4
         if      (indx2(i) .eq. inda) then 
            i2a = i 
         else if (indx2(i) .eq. indb) then
            i2b = i 
         else if (indx2(i) .eq. indc) then
            i2c = i
         else if (indx2(i) .eq. indd) then
            i2d = i
         else if (indx2(i) .eq. cind(1)) then
            i2p = i
         else if (indx2(i) .eq. cind(2)) then
            i2q = i
         else
            print *,'Error: Invalid index for x2 in cwork4442'
            print *,'X2 index is ',indx2(i),
     *              ' y indices: ',inda,indb,indc,indd,
     *              ' X2 indices: ',(indx2(j),j=1,4),
     *              ' cind = ',cind(1),cind(2)
            call abort_job() 
         end if
      end do

c   o make p and q canonical order
      if (i1p.gt.i1q) then
         i   = i1q
         i1q = i1p
         i1p = i
         i   = i2q
         i2q = i2p
         i2p = i
      end if
      b2pq = (i2p.lt.i2q)

      if      (i1p .eq. 1) then
         p1 = e1
         p2 = e2
      else if (i1p .eq. 2) then
         p1 = f1
         p2 = f2
      else if (i1p .eq. 3) then
         p1 = g1
         p2 = g2
      else if (i1p .eq. 4) then
         p1 = h1
         p2 = h2
      end if

      if      (i1q .eq. 1) then
         q1 = e1
         q2 = e2
      else if (i1q .eq. 2) then
         q1 = f1
         q2 = f2
      else if (i1q .eq. 3) then
         q1 = g1
         q2 = g2
      else if (i1q .eq. 4) then
         q1 = h1
         q2 = h2
      end if

      na = a2-a1+1
      nb = b2-b1+1
      nc = c2-c1+1
      nd = d2-d1+1
      ne = e2-e1+1
      nf = f2-f1+1
      ng = g2-g1+1
      nh = h2-h1+1
      ni = i2-i1+1
      nj = j2-j1+1
      nk = k2-k1+1
      nl = l2-l1+1
      np = p2-p1+1
      nq = q2-q1+1

      flopcount = 2*na*nb*nc*nd*np*nq

c------------------------------------------------------------------------
c   Design:
c      There are three main branches, each corresponds to X1 having the
c   b, c, or d output index. Within each branch, there are four sections:
c    1) sort X1
c    2) sort X2
c    3) dgemm
c    4) sort Y
c   Depending on which top-level branch is taken, scr1, scr2, and y will
c   be used differently.
c      The sorting sections are not exhaustive. As new SIAL code is
c   written, new permutations might have to be added. This is also by
c   design so a human can figure out the quickest way to sort the data.
c------------------------------------------------------------------------

      op_a = 'x'
      op_b = 'x'

c--------------------------------------------------------------------------
c B BRANCH : Y(a,b,c,d) = X1(a|b|(p,q)) * X2(c|d|p|q)
c--------------------------------------------------------------------------

      if (i1b.ne.0) then

c      o scr1: sorted X1 into either (a,b,p,q) or (p,q,a,b)
         if (i1a.eq.1) then
            if (i1b.eq.2) then
c            o x1(abpq) -> scr1(abpq)
               op_a = 'N'
               lda = na*nb
               n = na*nb*np*nq
               call xdcopy(n,x1,1,scr1,1)
            else if (i1b.eq.3) then
c-----------------------------------------------------------------------------
c   Transpose x1(a,p,b,q) --> scr1(a,b,p,q)
c-----------------------------------------------------------------------------

               n = 0
               do h = h1,h2
               do f = f1,f2
               do g = g1,g2
               do e = e1,e2
                  n = n + 1
                  scr1(n) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo
               
               op_a = 'N'
               lda = na*nb
c           else if (i1b.eq.4) then
            else
            end if
         else if (i1a.eq.2) then
            if (i1b.eq.1) then
c            o x1(b,a,pq) -> scr1(a,b,pq)
               op_a = 'N'
               lda = na*nb
               i = nb
               j = na
               n = np*nq
               call trnsps_bra(i,j,n,x1,i,j,scr1,j,i)
            else if (i1b.eq.3) then
c           else if (i1b.eq.4) then
            else
            end if
         else if (i1a.eq.3) then
            if (i1b.eq.1) then
            else if (i1b.eq.2) then
c           else if (i1b.eq.4) then
            else
c            o x1(pqab) -> scr1(pqab)
               op_a = 'T'
               lda = np*nq
               n = np*nq*na*nb
               call xdcopy(n,x1,1,scr1,1)
            end if
c        else if (i1a.eq.4) then
         else
            if (i1b.eq.1) then
c            o x1(b,pqa) -> scr1(pqa,b)
               op_a = 'T'
               lda = np*nq
               i = nb
               j = np*nq*na
               call trnsps(i,j,x1,i,scr1,j)
            else if (i1b.eq.2) then
c           else if (i1b.eq.3) then
            else
c            o x1(pq,b,a) -> scr1(pq,a,b)
               op_a = 'T'
               lda = np*nq
               i = nb
               j = na
               n = np*nq
               call trnsps_ket(i,j,n,x1,n,i,scr1,n,j)
            end if
         end if

c      o scr2: sorted X2 into either (c,d,p,q) or (p,q,c,d)
         if (i2c.eq.1) then
            if (i2d.eq.2) then
               if (b2pq) then
c               o x2(cdpq) -> scr2(cdpq)
                  op_b = 'T'
                  ldb = nc*nd
                  n = nc*nd*np*nq
                  call xdcopy(n,x2,1,scr2,1)
               else
c               o x2(cd,q,p) -> scr2(cd,p,q)
                  op_b = 'T'
                  ldb = nc*nd
                  i = nq
                  j = np
                  n = nc*nd
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               end if
            else if (i2d.eq.3) then
               if (b2pq) then
               else
               end if
c           else if (i2d.eq.4) then
            else
               if (b2pq) then
c               o x2(c,pq,d) -> scr2(c,d,pq)
                  op_b = 'T'
                  ldb = nc*nd
                  i = np*nq
                  j = nd
                  n = nc
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               else
c               o x2(c,qp,d) -> scr3(c,d,qp)
c               o scr3(cd,q,p) -> scr2(cd,p,q)
                  op_b = 'T'
                  ldb = nc*nd
                  i = nq*np
                  j = nd
                  n = nc
                  call trnsps_ket(i,j,n,x2,n,i,scr3,n,j)
                  i = nq
                  j = np
                  n = nc*nd
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               end if
            end if
         else if (i2c.eq.2) then
            if (i2d.eq.1) then
               if (b2pq) then
               else
               end if
            else if (i2d.eq.3) then
               if (b2pq) then
               else
               end if
c           else if (i2d.eq.4) then
            else
               if (b2pq) then
c----------------------------------------------------------------------------
c   Transpose x2(p,c,q,d) --> scr2(p,q,c,d)
c----------------------------------------------------------------------------
                  n = 0
                  do l = l1,l2
                  do j = j1,j2
                  do k = k1,k2
                  do i = i1,i2
                     n = n + 1
                     scr2(n) = x2(i,j,k,l)
                  enddo
                  enddo
                  enddo
                  enddo

                  op_b = 'N'
                  ldb = np*nq
               else
               end if
            end if
         else if (i2c.eq.3) then
            if (i2d.eq.1) then
               if (b2pq) then
               else
               end if
            else if (i2d.eq.2) then
               if (b2pq) then
c               o x2(p,dc,q) -> scr3(p,q,dc)
c               o scr3(pq,d,c) -> scr2(pq,c,d)
                  op_b = 'N'
                  ldb = np*nq
                  i = nd*nc
                  j = nq
                  n = np
                  call trnsps_ket(i,j,n,x2,n,i,scr3,n,j)
                  i = nd
                  j = nc
                  n = np*nq
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               else
c               o x2(qdc,p) -> scr3(p,qdc)
c               o scr3(pq,d,c) -> scr2(pq,c,d)
                  op_b = 'N'
                  ldb = np*nq
                  i = nq*nd*nc
                  j = np
                  call trnsps(i,j,x2,i,scr3,j)
                  i = nd
                  j = nc
                  n = np*nq
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               end if
c           else if (i2d.eq.4) then
            else
               if (b2pq) then
c               o x2(pqcd) -> scr2(pqcd)
                  op_b = 'N'
                  ldb = np*nq
                  n = np*nq*nc*nd
                  call xdcopy(n,x2,1,scr2,1)
               else
c               o x2(q,p,cd) -> scr2(p,q,cd)
                  op_b = 'N'
                  ldb = np*nq
                  i = nq
                  j = np
                  n = nc*nd
                  call trnsps_bra(i,j,n,x2,i,j,scr2,j,i)
               end if
            end if
c        else if (i2c.eq.4) then
         else
            if (i2d.eq.1) then
               if (b2pq) then
               else
               end if
            else if (i2d.eq.2) then
               if (b2pq) then
               else
               end if
c           else if (i2d.eq.3) then
            else
               if (b2pq) then
               else
               end if
            end if
         end if

         if (op_a.ne.'x'.and.op_b.ne.'x') then
c         o y(ab,cd) = scr1(ab|pq) * scr2(pq|cd)
            m = na*nb
            nn = nc*nd
            k = np*nq
            ldc = m
            call dgemm(op_a,op_b,m,nn,k,
     &                 1.d0,scr1,lda,
     &                      scr2,ldb,
     &                 0.d0,y,   ldc)
         else
            print *, '@CWORK4442: ADD MORE OPTIONS'
            print *,'         (i1a,i1b,op_a) = ',i1a,i1b,op_a
            print *,'     (i2c,i2d,i2p,op_b) = ',i2c,i2d,i2p,op_b
            call cwork4442_unopt(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     &                               a1, a2, b1, b2, c1, c2, d1, d2,
     &                           inda, indb, indc, indd,
     &                           x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     &                               e1, e2, f1, f2, g1, g2, h1, h2,
     &                           indx1,
     &                           x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     &                               i1, i2, j1, j2, k1, k2, l1, l2,
     &                           indx2,
     &                           cind, flopcount, scr1, scr2)
         end if

c--------------------------------------------------------------------------
c C BRANCH : Y(a,b,c,d) = X1(a|c|(p,q)) * X2(d|b|p|q)
c--------------------------------------------------------------------------

      else if (i1c.ne.0) then

c      o scr1: sorted X1 into either (a,c,p,q) or (p,q,a,c)
         if (i1a.eq.1) then
            if (i1c.eq.2) then
               op_a = 'N'
               lda = na*nc
               n = na*nc*np*nq
               call xdcopy(n,x1,1,scr1,1)
            else if (i1c.eq.3) then
c            o x1(a,p,c,q) -> scr1(a,c,p,q)
               op_a = 'N'
               lda = na*nc
               e = ne*ng*nf
               n = 1
               do h = h1,h2
                  call trnsps_ket(nf,ng,ne,x1(e1,f1,g1,h),ne,nf,
     &                                     scr1(n),ne,ng)
                  n = n + e
               end do
c           else if (i1c.eq.4) then
            else
c            o x1(a,pq,c) -> scr1(a,c,pq)
               op_a = 'N'
               lda = na*nc
               i = np*nq
               j = nc
               n = na
               call trnsps_ket(i,j,n,x1,n,i,scr1,n,j)
            end if
         else if (i1a.eq.2) then
            if (i1c.eq.1) then
            else if (i1c.eq.3) then
c           else if (i1c.eq.4) then
            else
c            o x1(p,a,q,c) -> scr1(p,q,a,c)
               op_a = 'T'
               lda = np*nq
               e = ne*ng*nf
               n = 1
               do h = h1,h2
                  call trnsps_ket(nf,ng,ne,x1(e1,f1,g1,h),ne,nf,
     &                                     scr1(n),ne,ng)
                  n = n + e
               end do
            end if
         else if (i1a.eq.3) then
            if (i1c.eq.1) then
c            o x1(cp,a,q) -> scr1(a,cp,q)
               op_a = 'N'
               lda = na*nc
               i = nc*np
               j = na
               n = nq
               call trnsps_bra(i,j,n,x1,i,j,scr1,j,i)
            else if (i1c.eq.2) then
c           else if (i1c.eq.4) then
            else
c            o x1(pqac) -> scr1(pqac)
               op_a = 'T'
               lda = np*nq
               n = np*nq*na*nc
               call xdcopy(n,x1,1,scr1,1)
            end if
c        else if (i1a.eq.4) then
         else
            if (i1c.eq.1) then
            else if (i1c.eq.2) then
c            o x1(p,c,qa) -> scr1(p,qa,c)
               op_a = 'T'
               lda = np*nq
               i = nc
               j = nq*na
               n = np
               call trnsps_ket(i,j,n,x1,n,i,scr1,n,j)
c           else if (i1c.eq.3) then
            else
            end if
         end if

c      o scr2: sorted X2 into either (d,b,p,q) or (p,q,d,b)
         if (i2b.eq.1) then
            if (i2d.eq.2) then
               if (b2pq) then
c               o x2(b,d,pq) -> scr2(d,b,pq)
                  op_b = 'T'
                  ldb = nd*nb
                  i = nb
                  j = nd
                  n = np*nq
                  call trnsps_bra(i,j,n,x2,i,j,scr2,j,i)
               else
               end if
            else if (i2d.eq.3) then
               if (b2pq) then
c               o x2(bp,d,q) -> scr2(d,bp,q)
                  op_b = 'T'
                  ldb = nd*nb
                  i = nb*np
                  j = nd
                  n = nq
                  call trnsps_bra(i,j,n,x2,i,j,scr2,j,i)
               else
c               o x2(b,q,dp) -> scr3(b,dp,q)
c               o scr3(b,d,pq) -> scr2(d,b,pq)
                  op_b = 'T'
                  ldb = nd*nb
                  i = nq
                  j = nd*np
                  n = nb
                  call trnsps_ket(i,j,n,x2,n,i,scr3,n,j)
                  i = nb
                  j = nd
                  n = np*nq
                  call trnsps_bra(i,j,n,scr3,i,j,scr2,j,i)
               end if
c           else if (i2d.eq.4) then
            else
               if (b2pq) then
c               o x2(b,pqd) -> scr2(pqd,b)
                  op_b = 'N'
                  ldb = np*nq
                  i = nb
                  j = np*nq*nd
                  call trnsps(i,j,x2,i,scr2,j)
               else
               end if
            end if
         else if (i2b.eq.2) then
            if (i2d.eq.1) then
               if (b2pq) then
               else
               end if
            else if (i2d.eq.3) then
               if (b2pq) then
               else
               end if
c           else if (i2d.eq.4) then
            else
               if (b2pq) then
c               o x2(p,b,qd) -> scr2(p,qd,b)
                  op_b = 'N'
                  ldb = np*nq
                  i = nb
                  j = nq*nd
                  n = np
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               else
c               o x2(qb,p,d) -> scr3(p,qb,d)
c               o scr3(pq,b,d) -> scr2(pq,d,b)
                  op_b = 'N'
                  ldb = np*nq
                  i = nq*nb
                  j = np
                  n = nd
                  call trnsps_bra(i,j,n,x2,i,j,scr3,j,i)
                  i = nb
                  j = nd
                  n = np*nq
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               end if
            end if
         else if (i2b.eq.3) then
            if (i2d.eq.1) then
               if (b2pq) then
c               o x2(d,p,b,q) -> scr2(d,b,p,q)
                  op_b = 'T'
                  ldb = nd*nb
                  i = ni*nk*nj
                  n = 1
                  do l = l1,l2
                     call trnsps_ket(nj,nk,ni,x2(i1,j1,k1,l),ni,nj,
     &                                        scr2(n),ni,nk)
                     n = n + i
                  end do
               else
c               o x2(d,q,bp) -> scr2(d,bp,q)
                  op_b = 'T'
                  ldb = nd*nb
                  i = nq
                  j = nb*np
                  n = nd
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               end if
            else if (i2d.eq.2) then
               if (b2pq) then
               else
               end if
c           else if (i2d.eq.4) then
            else
               if (b2pq) then
c               o x2(pq,b,d) -> scr2(pq,d,b)
                  op_b = 'N'
                  ldb = np*nq
                  i = nb
                  j = nd
                  n = np*nq
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               else
               end if
            end if
c        else if (i2b.eq.4) then
         else
            if (i2d.eq.1) then
               if (b2pq) then
               else
               end if
            else if (i2d.eq.2) then
               if (b2pq) then
c               o x2(p,d,q,b) -> scr2(p,q,d,b)
                  op_b = 'N'
                  ldb = np*nq
                  i = ni*nk*nj
                  n = 1
                  do l = l1,l2
                     call trnsps_ket(nj,nk,ni,x2(i1,j1,k1,l),ni,nj,
     &                                        scr2(n),ni,nk)
                     n = n + i
                  end do
               else
               end if
c           else if (i2d.eq.3) then
            else
               if (b2pq) then
c              o x2(p,q,d,b) -> scr2(p,q,d,b)
                  op_b = 'N'
                  ldb = np*nq
                  n = np*nq*nb*nd
                  call xdcopy(n,x2,1,scr2,1)
               else
c              o x2(q,p,d,b) -> scr2(p,q,d,b)
               end if
            end if
         end if

         if (op_a.ne.'x'.and.op_b.ne.'x') then
c         o scr3(ad,bc) = scr1(ac|pq) * scr2(pq|db)
            m = na*nc
            nn = nd*nb
            k = np*nq
            ldc = m
            call dgemm(op_a,op_b,m,nn,k,
     &                 1.d0,scr1,lda,
     &                      scr2,ldb,
     &                 0.d0,scr3,ldc)
c         o scr3(a,cd,b) -> y(a,b,cd)
            i = nc*nd
            j = nb
            n = na
            call trnsps_ket(i,j,n,scr3,n,i,y,n,j)
         else
            print *,'@CWORK4442: ADD MORE OPTIONS'
            print *,'          (i1a,i1c,op_a) = ',i1a,i1c,op_a
            print *,'      (i2b,i2d,i2p,op_b) = ',i2b,i2d,i2p,op_b
            call cwork4442_unopt(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     &                               a1, a2, b1, b2, c1, c2, d1, d2,
     &                           inda, indb, indc, indd,
     &                           x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     &                               e1, e2, f1, f2, g1, g2, h1, h2,
     &                           indx1,
     &                           x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     &                               i1, i2, j1, j2, k1, k2, l1, l2,
     &                           indx2,
     &                           cind, flopcount, scr1, scr2)
         end if

c--------------------------------------------------------------------------
c D BRANCH : Y(a,b,c,d) = X1(a|d|(p,q)) * X2(b|c|p|q)
c--------------------------------------------------------------------------

      else if (i1d.ne.0) then

c      o scr1: sorted X1 into either (a,d,p,q) or (p,q,a,d)
         if (i1a.eq.1) then
            if (i1d.eq.2) then
            else if (i1d.eq.3) then
c           else if (i1d.eq.4) then
            else
c            o x1(a,pq,d) -> scr1(a,d,pq)
               op_a = 'N'
               lda = na*nd
               i = np*nq
               j = nd
               n = na
               call trnsps_ket(i,j,n,x1,n,i,scr1,n,j)
            end if
         else if (i1a.eq.2) then
            if (i1d.eq.1) then
c---------------------------------------------------------------------------
c   Transpose x1(d,a,p,q) --> scr1(a,d,p,q)
c---------------------------------------------------------------------------
               op_a = 'N'
               lda  = na*nd
               i = nd
               j = na
               n = np*nq
               call trnsps_bra(i,j,n,x1,i,j,scr1,j,i)
            else if (i1d.eq.3) then
c            o x1(p,ad,q) -> scr1(p,q,ad)
               op_a = 'T'
               lda = np*nq
               i = na*nd
               j = nq
               n = np
               call trnsps_ket(i,j,n,x1,n,i,scr1,n,j)
c           else if (i1d.eq.4) then
            else
            end if
         else if (i1a.eq.3) then
            if (i1d.eq.1) then
            else if (i1d.eq.2) then
c            o x1(p,da,q) -> scr2(p,q,da)
c            o scr2(pq,d,a) -> scr1(pq,a,d)
               op_a = 'T'
               lda = np*nq
               i = nd*na
               j = nq
               n = np
               call trnsps_ket(i,j,n,x1,n,i,scr2,n,j)
               i = nd
               j = na
               n = np*nq
               call trnsps_ket(i,j,n,scr2,n,i,scr1,n,j)
c           else if (i1d.eq.4) then
            else
c            o x1(pqad) -> scr1(pqad)
               op_a = 'T'
               lda = np*nq
               n = np*nq*na*nd
               call xdcopy(n,x1,1,scr1,1)
            end if
c        else if (i1a.eq.4) then
         else
            if (i1d.eq.1) then
c            o x1(d,pqa) -> scr1(pqa,d)
               op_a = 'T'
               lda = np*nq
               i = nd
               j = np*nq*na
               call trnsps(i,j,x1,i,scr1,j)
            else if (i1d.eq.2) then
c           else if (i1d.eq.3) then
            else
c            o x1(pq,d,a) -> scr1(pq,a,d)
               op_a = 'T'
               lda = np*nq
               i = nd
               j = na
               n = np*nq
               call trnsps_ket(i,j,n,x1,n,i,scr1,n,j)
            end if
         end if

c      o scr2: sorted X2 into either (b,c,p,q) or (p,q,b,c)
         if (i2b.eq.1) then
            if (i2c.eq.2) then
               if (b2pq) then
c               o x2(bcpq) -> scr2(bcpq)
                  op_b = 'T'
                  ldb = nb*nc
                  n = nb*nc*np*nq
                  call xdcopy(n,x2,1,scr2,1)
               else
               end if
            else if (i2c.eq.3) then
               if (b2pq) then
               else
               end if
c           else if (i2c.eq.4) then
            else
               if (b2pq) then
c-----------------------------------------------------------------------------
c   Transpose x2(b,p,q,c) --> scr2(b,c,p,q)
c-----------------------------------------------------------------------------

                  op_b = 'T'
                  ldb = nb*nc
                  i   = np*nq
                  j   = nc
                  n   = nb
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               else
               end if
            end if
         else if (i2b.eq.2) then
            if (i2c.eq.1) then
               if (b2pq) then
               else
c---------------------------------------------------------------------------
c   Transpose x2(c,b,q,p) --> scr2(b,c,p,q)
c---------------------------------------------------------------------------

                  n = 0
                  do k = k1,k2
                  do l = l1,l2
                  do i = i1,i2
                  do j = j1,j2
                     n = n + 1
                     scr2(n) = x2(i,j,k,l)
                  enddo
                  enddo
                  enddo
                  enddo

                  op_b = 'T'
                  ldb = nb*nc 
               end if
            else if (i2c.eq.3) then
               if (b2pq) then
c               o x2(p,bc,q) -> scr2(p,q,bc)
                  op_b = 'N'
                  ldb = np*nq
                  i = nb*nc
                  j = nq
                  n = np
                  call trnsps_ket(i,j,n,x2,n,i,scr2,n,j)
               else
c               o x2(q,bcp) -> scr2(bcp,q)
                  op_b = 'T'
                  ldb = nb*nc
                  i = nq
                  j = nb*nc*np
                  call trnsps(i,j,x2,i,scr2,j)
               end if
c           else if (i2c.eq.4) then
            else
               if (b2pq) then
               else
               end if
            end if
         else if (i2b.eq.3) then
            if (i2c.eq.1) then
               if (b2pq) then
               else
               end if
            else if (i2c.eq.2) then
               if (b2pq) then
c               o x2(p,cb,q) -> scr3(p,q,cb)
c               o scr3(pq,c,b) -> scr2(pq,b,c)
                  op_b = 'N'
                  ldb = np*nq
                  i = nc*nb
                  j = nq
                  n = np
                  call trnsps_ket(i,j,n,x2,n,i,scr3,n,j)
                  i = nc
                  j = nb
                  n = np*nq
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               else
c               o x2(qcb,p) -> scr3(p,qcb)
c               o scr3(pq,c,b) -> scr2(pq,b,c)
                  op_b = 'N'
                  ldb = np*nq
                  i = nq*nc*nb
                  j = np
                  call trnsps(i,j,x2,i,scr3,j)
                  i = nc
                  j = nb
                  n = np*nq
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               end if
c           else if (i2c.eq.4) then
            else
               if (b2pq) then
c               o x2(pqbc) -> scr2(pqbc)
                  op_b = 'N'
                  ldb = np*nq
                  n = np*nq*nb*nc
                  call xdcopy(n,x2,1,scr2,1)
               else
               end if
            end if
c        else if (i2b.eq.4) then
         else
            if (i2c.eq.1) then
               if (b2pq) then
c               o x2(c,pqb) -> scr2(pqb,c)
                  op_b = 'N'
                  ldb = np*nq
                  i = nc
                  j = np*nq*nb
                  call trnsps(i,j,x2,i,scr2,j)
               else
c               o x2(cqp,b) -> scr3(b,cqp)
c               o scr3(bc,q,p) -> scr2(bc,p,q)
                  op_b = 'T'
                  ldb = nb*nc
                  i = nc*nq*np
                  j = nb
                  call trnsps(i,j,x2,i,scr3,j)
                  i = nq
                  j = np
                  n = nb*nc
                  call trnsps_ket(i,j,n,scr3,n,i,scr2,n,j)
               end if
            else if (i2c.eq.2) then
               if (b2pq) then
               else
               end if
c           else if (i2c.eq.3) then
            else
               if (b2pq) then
               else
               end if
            end if
         end if

         if (op_a.ne.'x'.and.op_b.ne.'x') then
c         o scr3(ad,bc) = scr1(ad|pq) * scr2(pq|bc)
            m = na*nd
            nn = nb*nc
            k = np*nq
            ldc = m
            call dgemm(op_a,op_b,m,nn,k,
     &                 1.d0,scr1,lda,
     &                      scr2,ldb,
     &                 0.d0,scr3,ldc)
c         o scr3(a,d,bc) -> y(a,bc,d)
            i = nd
            j = nb*nc
            n = na
            call trnsps_ket(i,j,n,scr3,n,i,y,n,j)
         else
            print *,'@CWORK4442: ADD MORE OPTIONS'
            print *,'          (i1a,i1d,op_a) = ',i1a,i1d,op_a
            print *,'      (i2b,i2c,i2p,op_b) = ',i2b,i2c,i2p,op_b
            call cwork4442_unopt(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     &                               a1, a2, b1, b2, c1, c2, d1, d2,
     &                           inda, indb, indc, indd,
     &                           x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     &                               e1, e2, f1, f2, g1, g2, h1, h2,
     &                           indx1,
     &                           x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     &                               i1, i2, j1, j2, k1, k2, l1, l2,
     &                           indx2,
     &                           cind, flopcount, scr1, scr2)
         end if

c     end if (three top-level branches)
      end if

c ----------------------------------------------------------------------

      if (debug) then
         if (op_a.ne.'x'.and.op_b.ne.'x') then

c-----------------------------------------------------------------------
c   Save the results of the optimized contraction.
c-----------------------------------------------------------------------

            n = na*nb*nc*nd
            call xdcopy(n, y, 1, scr1, 1)

c-----------------------------------------------------------------------
c   Perform the contraction again in unoptimized mode.
c-----------------------------------------------------------------------

            call cwork4442_unopt(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     &                               a1, a2, b1, b2, c1, c2, d1, d2,
     &                           inda, indb, indc, indd,
     &                           x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     &                               e1, e2, f1, f2, g1, g2, h1, h2,
     &                           indx1,
     &                           x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     &                               i1, i2, j1, j2, k1, k2, l1, l2,
     &                           indx2,
     &                           cind, flopcount, scr1, scr2)

c-----------------------------------------------------------------------
c   Compare the results.
c-----------------------------------------------------------------------

            do i = 1, n
               if (dabs(y(i)-scr1(i)) .gt. 1.d-12) then
                  print *,'Error in cwork4442: contractions do not ',
     *                      'match at line ',current_line
                  print *,'OPT: ',scr1(i),' UNOPT: ',y(i)
                  print *,'a1,a2,b1,b2,c1,c2,d1,d2 ',
     *                     a1,a2,b1,b2,c1,c2,d1,d2
                  print *,'e1,e2,f1,f2,g1,g2,h1,h2 ',
     *                     e1,e2,f1,f2,g1,g2,h1,h2
                  print *,'i1,i2,j1,j2,k1,k2,l1,l2 ',
     *                     i1,i2,j1,j2,k1,k2,l1,l2
                  print *,'indx1: ',(indx1(a),a=1,4)
                  print *,'indx2: ',(indx2(a),a=1,4)
                  call abort_job()
               end if 
            end do

c        end if (optimized cont)
         end if
c     end if (debug)
      end if

      return
      end


