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
       subroutine cwork4643(y, na1,na2,nb1,nb2,
     *                         nc1,nc2,nd1,nd2,
     *                         a1,a2,b1,b2,
     *                         c1, c2, d1, d2,
     *                         inda, indb, indc, indd,
c
     *                      x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *                         ne1,ne2,nf1,nf2,
     *                          m1, m2, n1, n2, r1, r2, s1, s2,
     *                          e1, e2, f1, f2,
     *                      indx1,
c
     *                      x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                          i1, i2, j1, j2, k1, k2, l1, l2,
     *                      indx2,
     *                      cind, flopcount, scr1, scr2, scr3)
c-------------------------------------------------------------------------
c   Performs a "4643" contraction: 
c      4 index output array  : y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
c      6 index operand array : x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,ne1:ne2,nf1:nf2)
c      4 index operand array : x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
c      3 index contraction.
c
c   I. e., contract two of the the 4 indices of the 2nd operand array out 
c   of the first operand array, replacing the indices by the 2nd operand's 
c   non-contracted indices.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *        ne1,ne2,nf1,nf2,
     *        nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2, 
     *        ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2
      integer a1,a2,b1,b2,c1,c2,d1,d2,
     *        e1,e2,f1,f2,
     *        n1,n2,m1,m2,r1,r2,s1,s2,
     *        i1,i2,j1,j2,k1,k2,l1,l2
      integer inda, indb, indc, indd, inde, indf, indx1(6), indx2(4)
      integer cind(3), flopcount

      double precision y(*)
      double precision x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,
     *                    ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
      double precision scr1(*), scr2(*), scr3(*)
      double precision etemp  
c
c Temporary arrays used in array permutations 
c 
c END Temporary arrays used in array permutations 

c     integer p1, p2, q1, q2
      integer na, nb, nc, nd, ne, nf, ng, nh, ni, nj, nk, nl, np, nq
      integer p0(4), q1(6), q2(4) 

      integer i1a, i1b, i1c, i1d, i1p, i1q
      integer i2a, i2b, i2c, i2d, i2p, i2q
      integer a, b, c, d, e, f, g, h, i, j, l, p, q
      integer r, s 

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
      if (m1.ne.nm1.or.n1.ne.nn1.or.r1.ne.nr1.or.s1.ne.ns1.or.
     &    m2.ne.nm2.or.n2.ne.nn2.or.r2.ne.nr2.or.s2.ne.ns2.or.    
     &    e1.ne.ne1.or.e2.ne.ne2.or.f1.ne.nf1.or.f2.ne.nf2    ) then
         print *, '@CWORK4643: Assertion failed.'
         print *, '            The x1 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert x2 is contiguous
      if (i1.ne.ni1.or.j1.ne.nj1.or.k1.ne.nk1.or.l1.ne.nl1.or.
     &    i2.ne.ni2.or.j2.ne.nj2.or.k2.ne.nk2.or.l2.ne.nl2 ) then 
         print *, '@CWORK4643: Assertion failed.'
         print *, '            The x2 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert y is contiguous
      if (a1.ne.na1.or.b1.ne.nb1.or.c1.ne.nc1.or.d1.ne.nd1) then 
         print *, '@CWORK4643: Assertion failed.'
         print *, '            The product tensor is not contiguous.'
         call abort_job()
      end if

c     end if (.false.)
      end if

c-----------------------------------------------------------------------
c   Perform the contraction in unoptimized mode.
c-----------------------------------------------------------------------
            call cwork4643_unopt(y, na1,na2,nb1,nb2, 
     *                         nc1,nc2,nd1,nd2, 
     *                         a1,a2,b1,b2,  
     *                         c1, c2, d1, d2,  
     *                         inda, indb, indc, indd,   
c 
     *                      x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *                         ne1,ne2,nf1,nf2, 
     *                          m1, m2, n1, n2, r1, r2, s1, s2, 
     *                          e1, e2, f1, f2, 
     *                      indx1,
c 
     *                      x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                          i1, i2, j1, j2, k1, k2, l1, l2,
     *                      indx2,
     *                      cind, flopcount, scr1, scr2, scr3) 

      return
      end


