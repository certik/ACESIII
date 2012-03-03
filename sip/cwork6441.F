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
       subroutine cwork6441(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                         ne1,ne2,nf1,nf2,
     *                      a1, a2, b1, b2, c1, c2, d1, d2,
     *                      e1, e2, f1, f2,
     *                      inda, indb, indc, indd, inde, indf,
c
     *                      x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *                          m1, m2, n1, n2, r1, r2, s1, s2,
     *                      indx1,
c
     *                      x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                          i1, i2, j1, j2, k1, k2, l1, l2,
     *                      indx2,
     *                      cind, flopcount, scr1, scr2, scr3)
c      6 index output array  : y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,ne1:ne2,nf1:nf2)
c      4 index operand array : x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2)
c      4 index operand array : x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
c      1 index contraction.
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
      integer inda, indb, indc, indd, inde, indf, indx1(4), indx2(4)
      integer cind(1), flopcount

c     double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,
c    *                   ne1:ne2,nf1:nf2)
      double precision y(*)
      double precision x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
      double precision scr1(*), scr2(*), scr3(*)
      double precision etemp, thresh, max1, max2   
c
c Temporary arrays used in array permutations 
c 
c
c END Temporary arrays used in array permutations 

c     integer p1, p2, q1, q2
      integer na, nb, nc, nd, ne, nf, ng, nh, ni, nj, nk, nl, np, nq
      integer p0(6), q1(4), q2(4) 

      integer i1a, i1b, i1c, i1d, i1p, i1q
      integer i2a, i2b, i2c, i2d, i2p, i2q
      integer a, b, c, d, e, f, g, h, i, j, l, p, q
      integer aa, bb, cc, dd, ee, ff, ii, jj, ll  
      integer yy, xx, y1, y2, y3, y4, yTOT, nTOT, inc, xinc    
      integer ystart, yend 

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
     &    m2.ne.nm2.or.n2.ne.nn2.or.r2.ne.nr2.or.s2.ne.ns2    ) then
         print *, '@CWORK6441: Assertion failed.'
         print *, '            The x1 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert x2 is contiguous
      if (i1.ne.ni1.or.j1.ne.nj1.or.k1.ne.nk1.or.l1.ne.nl1.or.
     &    i2.ne.ni2.or.j2.ne.nj2.or.k2.ne.nk2.or.l2.ne.nl2    ) then
         print *, '@CWORK6441: Assertion failed.'
         print *, '            The x2 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert y is contiguous
      if (a1.ne.na1.or.b1.ne.nb1.or.c1.ne.nc1.or.d1.ne.nd1.or.
     &    a2.ne.na2.or.b2.ne.nb2.or.c2.ne.nc2.or.d2.ne.nd2.or. 
     &    e1.ne.ne1.or.f1.ne.nf1.or.e2.ne.ne2.or.f2.ne.nf2    ) then
         print *, '@CWORK6441: Assertion failed.'
         print *, '            The product tensor is not contiguous.'
         call abort_job()
      end if

c     end if (.false.)
      end if

c     go to 99 
c Try optimized version !!! VFL 

      inc = 5 

      do i = 1, 6 
         p0(i) = i 
      enddo 
      do i = 1, 4 
         q1(i) = 0 
         q2(i) = 0 
      enddo 

      do i = 1, 4 
         if (indx1(i) .eq. inda) q1(i) = p0(1) 
         if (indx1(i) .eq. indb) q1(i) = p0(2) 
         if (indx1(i) .eq. indc) q1(i) = p0(3) 
         if (indx1(i) .eq. indd) q1(i) = p0(4) 
         if (indx1(i) .eq. inde) q1(i) = p0(5) 
         if (indx1(i) .eq. indf) q1(i) = p0(6) 

         if (indx2(i) .eq. inda) q2(i) = p0(1) 
         if (indx2(i) .eq. indb) q2(i) = p0(2) 
         if (indx2(i) .eq. indc) q2(i) = p0(3) 
         if (indx2(i) .eq. indd) q2(i) = p0(4) 
         if (indx2(i) .eq. inde) q2(i) = p0(5) 
         if (indx2(i) .eq. indf) q2(i) = p0(6) 
      enddo 

      do i = 1, 4 
      do j = 1, 4 
         if (indx1(i) .eq. indx2(j)) then 
            q1(i) = 7 ! indx1(i) 
            q2(j) = 7 ! indx1(i) 
         endif 
      enddo 
      enddo 

      nTOT = (f2-f1+1) 
     *     * (e2-e1+1) 
     *     * (d2-d1+1) 
     *     * (c2-c1+1) 
     *     * (b2-b1+1) 
     *     * (a2-a1+1) 

      thresh = 100.0 ! 1.0d-10 

      do a  = 1, nTOT   
         y(a) = 0.0 
      enddo 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.7).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.7).and.(q2(3).eq.5).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
            y2 = 0 
            do c = c1, c2 
            do j = j1, j2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,j,e,f) 
            enddo 
            enddo 

         do d = d1, d2
            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do j = j1, j2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,j,d) 
            enddo 
            enddo 
            enddo 
            y4 = 0   
         do c = c1, c2  
            y4 = y4 + 1 
            xx = (j2-j1+1)*(y4-1) 
            y1 = 0 

         do b = b1, b2, inc 
            do bb = b, min(b2,b+inc-1) 
         do a = a1, a2, inc 
            do aa = a, min(a2,a+inc-1) 

               etemp = 0.0 
               yy    = yy + 1 

               y2 = xx  
               do j = j1, j2     
                  y1 = y1 + 1 
                  y2 = y2 + 1   
                  etemp = etemp + scr1(y1)*scr2(y2)
               enddo

               y(yy) = etemp 
         
            enddo 
         enddo
            enddo 
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.6).and.(q1(3).eq.7).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.7).and.(q2(3).eq.5).and.
     *    (q2(4).eq.2)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do a = a1, a2  
            do j = j1, j2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,f,j,d) 
            enddo 
            enddo 
         do c = c1, c2  
         do b = b1, b2  
            y1 = 0 
            y2 = 0 
            do j = j1, j2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,j,e,b) 
            enddo 
         do a = a1, a2  

            etemp = 0.0 

            y2 = 0 
            do j = j1, j2 
               y1 = y1 + 1 
               y2 = y2 + 1 
c              etemp = etemp + x1(a,f,j,d)*x2(c,j,e,b)
               etemp = etemp + scr1(y1)*scr2(y2)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.5).and.(q1(2).eq.2).and.(q1(3).eq.7).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.7).and.(q2(3).eq.1).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0.0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2

            y2 = 0 
            do a = a1, a2
            do j = j1, j2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,j,a,f) 
            enddo 
            enddo 

         do b = b1, b2

            y2 = 0 
            y1 = 0 
            do j = j1, j2 
               y1 = y1 + 1 
               scr1(y1) = x1(e,b,j,d) 
            enddo 

         do a = a1, a2

            etemp = 0.0 

            y1 = 0 
            do j = j1, j2 
c              etemp = etemp + x1(e,b,j,d)*x2(c,j,a,f)
               y2 = y2 + 1 
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(y2)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.7).and.(q1(3).eq.5).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
            y1 = 0  
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,k,e,f) 
            enddo 
            enddo 
         do d = d1, d2
         do c = c1, c2
         do b = b1, b2
            y1 = 0 
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,b,k,d) 
            enddo 
         do a = a1, a2, inc 
            do aa = a, min(a2,a+inc-1) 

               etemp = 0.0 

               y2 = 0 
               do k = k1, k2 
                  y1 = y1 + 1 
                  y2 = y2 + 1 
c                 etemp = etemp + x1(a,k,e,f)*x2(c,b,k,d)
                  etemp = etemp + scr1(y1)*scr2(y2)
               enddo

               yy = yy + 1 
               y(yy) = etemp  
         
         enddo
            enddo 
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.7).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.7).and.(q2(3).eq.5).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2

            y1 = 0 
            do b = b1, b2 
            do a = a1, a2 
            do j = j1, j2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,j,f) 
            enddo 
            enddo 
            enddo 

         do e = e1, e2
         do d = d1, d2
         do c = c1, c2

            y2 = 0 
            do j = j1, j2 
               y2 = y2 + 1  
               scr2(y2) = x2(c,j,e,d) 
            enddo 

            y1 = 0 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 
 
            do j = j1, j2 
c              etemp = etemp + x1(a,b,j,f)*x2(c,j,e,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(j-j1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.7).and.(q1(3).eq.5).and.
     *    (q1(4).eq.2).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.6).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2

            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,k,e,b) 
            enddo 
            enddo 
            enddo 

         do d = d1, d2
         do c = c1, c2
            y1 = 0 

            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,f,k,d) 
            enddo 

         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,k,e,b)*x2(c,f,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.4).and.(q1(3).eq.5).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2

            y1 = 0  
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,d,e,k) 
            enddo 
            enddo 

         do c = c1, c2
         do b = b1, b2
            y1 = 0 

            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,b,k,f) 
            enddo 

         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,d,e,k)*x2(c,b,k,f)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.4).and.(q1(3).eq.3).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.7).and.(q2(2).eq.2).and.(q2(3).eq.5).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2

            y2 = 0 
            do b = b1, b2
            do i = i1, i2 
               y2 = y2 + 1 
               scr2(y2) = x2(i,b,e,f) 
            enddo 
            enddo 

         do d = d1, d2
         do c = c1, c2

            y1 = 0 
            do a = a1, a2
            do i = i1, i2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,d,c,i) 
            enddo 
            enddo 

            y4 = 0 
         do b = b1, b2
            y4 = y4 + 1 
            xx = (i2-i1+1)*(y4-1) 
            y1 = 0 
         do a = a1, a2

               etemp = 0.0 
               y2 = xx  
            do i = i1, i2 
c              etemp = etemp + x1(a,d,c,i)*x2(i,b,e,f)
               y1 = y1 + 1 
               y2 = y2 + 1 
               etemp = etemp + scr1(y1)*scr2(y2)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.7).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.4).and.(q2(3).eq.5).and.
     *    (q2(4).eq.7)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2

            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do l = l1, l2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,l,f) 
            enddo 
            enddo 
            enddo 

         do e = e1, e2
         do d = d1, d2
         do c = c1, c2

            y2 = 0 
            do l = l1, l2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,d,e,l) 
            enddo 

            y1 = 0 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do l = l1, l2 
c              etemp = etemp + x1(a,b,l,f)*x2(c,d,e,l)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(l-l1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.5).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.4).and.(q2(3).eq.7).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0
         do f = f1, f2
         do e = e1, e2

            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,e,k) 
            enddo 
            enddo 
            enddo 

         do d = d1, d2
         do c = c1, c2

            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,d,k,f) 
            enddo 

         y1 = 0 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,b,e,k)*x2(c,d,k,f)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.3).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.7).and.(q2(2).eq.2).and.(q2(3).eq.5).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
            y2 = 0 
            do b = b1, b2 
            do i = i1, i2 
               y2 = y2 + 1 
               scr2(y2) = x2(i,b,e,f) 
            enddo 
            enddo 

         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            do a = a1, a2 
            do i = i1, i2 
               y1 = y1 + 1 
               scr1(y1) = x1(c,i,a,d) 
            enddo 
            enddo 

         do b = b1, b2
            y1 = 0 
         do a = a1, a2

            etemp = 0.0 

            y3 = 0 
            do i = i1, i2 
               y1 = y1 + 1 
               y3 = y3 + 1 
               y2 = y3 + (b-b1)*(i2-i1+1)  
c              etemp = etemp + x1(c,i,a,d)*x2(i,b,e,f)
               etemp = etemp + scr1(y1)*scr2(y2)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo ! a 
         enddo ! b 
         enddo ! c 
         enddo ! d 
         enddo ! e 
         enddo ! f 

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.3).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.2).and.  
     *    (q2(1).eq.7).and.(q2(2).eq.4).and.(q2(3).eq.5).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
         do i = i1, i2 
            y1       = y1 + 1 
            scr1(y1) = x1(c,i,a,b) 
         enddo 
         enddo 
         enddo 
         enddo 

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y2 = 0 
            do i = i1, i2 
               y2       = y2 + 1 
               scr2(y2) = x2(i,d,e,f) 
            enddo 
            y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            y2 = 0 
            do i = i1, i2 
               y2 = y2 + 1 
               y1 = y1 + 1 
c              etemp = etemp + x1(c,i,a,b)*x2(i,d,e,f)
               etemp = etemp + scr1(y1)*scr2(y2)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.5).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(e,k,a,d) 
            enddo 
            enddo 

         do c = c1, c2
         do b = b1, b2
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,b,k,f) 
            enddo 
            y1 = 0 
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
               y1 = y1 + 1 
c              etemp = etemp + x1(e,k,a,d)*x2(c,b,k,f)
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1  
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.5).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.2).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.4).and.(q2(3).eq.7).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
            y2 = 0 
            do d = d1, d2
            do c = c1, c2
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,d,k,f) 
            enddo 
            enddo 
            enddo 
         do e = e1, e2
            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do k = k1, k2 
               y1       = y1 + 1 
               scr1(y1) = x1(e,k,a,b) 
            enddo 
            enddo 
            enddo 
         y4 = 0 
         do d = d1, d2
         do c = c1, c2
            y4 = y4 + 1 
            y3 = (y4-1) *(k2-k1+1)  
            y1 = 0 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            y2 = 0 
            do k = k1, k2 
               y1 = y1 + 1 
               y2 = y2 + 1 
c              etemp = etemp + x1(e,k,a,b)*x2(c,d,k,f)
               etemp = etemp + scr1(y1)*scr2(y2+y3)
            enddo

            yy = yy + 1  
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.4).and.(q1(3).eq.7).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.7).and.(q2(3).eq.3).and.
     *    (q2(4).eq.2)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2  
         do e = e1, e2, inc 
            do ee = e, min(e2,e+inc-1) 
         do d = d1, d2  
         do c = c1, c2  
         do b = b1, b2  
            y2 = 0 
            do j = j1, j2 
               y2 = y2 + 1 
               scr2(y2) = x2(ee,j,c,b) 
            enddo 
         do a = a1, a2,inc 
            do aa = a, min(a2,a+inc-1) 

               etemp = 0.0 

               do j = j1, j2   
c                 etemp = etemp + x1(aa,d,j,f)*x2(ee,j,c,b)
                  etemp = etemp + x1(aa,d,j,f)*scr2(j-j1+1)
               enddo

               yy = yy + 1 
               y(yy) = etemp  
         
            enddo 
         enddo
         enddo
         enddo
         enddo
            enddo 
         enddo
         enddo

         return 
c        go to 99 

      endif 

c VFL 
      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.7).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.7).and.(q2(3).eq.3).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         y3 = 0 
         do f = f1, f2
            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do j = j1, j2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,j,f) 
            enddo 
            enddo 
            enddo 
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
            y2 = 0 
            do j = j1, j2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,j,c,d) 
            enddo 
         y1 = 0 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            y2 = 0 
            do j = j1, j2 
               y1 = y1 + 1 
               y2 = y2 + 1 
c              etemp = etemp + x1(a,b,j,f)*x2(e,j,c,d)
               etemp = etemp + scr1(y1)*scr2(y2)  
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.4).and.(q1(3).eq.7).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.2).and.(q2(3).eq.5).and.
     *    (q2(4).eq.7)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do a = a1, a2
            do l = l1, l2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,d,l,f) 
            enddo 
            enddo 

         do c = c1, c2
         do b = b1, b2
            y1 = 0 

            y2 = 0 
            do l = l1, l2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,b,e,l) 
            enddo 

         do a = a1, a2

            etemp = 0.0 

            do l = l1, l2 
c              etemp = etemp + x1(a,d,l,f)*x2(c,b,e,l)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(l-l1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.3).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.7).and.(q2(2).eq.4).and.(q2(3).eq.5).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
         do i = i1, i2 
            y1 = y1 + 1 
            scr1(y1) = x1(a,b,c,i) 
         enddo 
         enddo 
         enddo 
         enddo 

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y2 = 0 
            do i = i1, i2 
               y2 = y2 + 1 
               scr2(y2) = x2(i,d,e,f) 
            enddo 
            y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do i = i1, i2 
c              etemp = etemp + x1(a,b,c,i)*x2(i,d,e,f)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(i-i1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.6).and.(q1(3).eq.5).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2

            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,f,e,k) 
            enddo 
            enddo 

         do d = d1, d2
         do c = c1, c2
         do b = b1, b2
            y1 = 0 

            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,b,k,d) 
            enddo 

         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,f,e,k)*x2(c,b,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 
 

      if ((q1(1).eq.1).and.(q1(2).eq.4).and.(q1(3).eq.5).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.6).and.(q2(3).eq.7).and.
     *    (q2(4).eq.2)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,d,e,k) 
            enddo 
            enddo 
         do c = c1, c2
         do b = b1, b2
            y1 = 0 
            y2 = 0  
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,f,k,b) 
            enddo 
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
               etemp = etemp + x1(a,d,e,k)*x2(c,f,k,b)
               y1 = y1 + 1 
c              etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.7).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.6).and.(q2(3).eq.3).and.
     *    (q2(4).eq.7)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do l = l1, l2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,l,d) 
            enddo 
            enddo 
            enddo 
           
         do c = c1, c2
            y1 = 0 
            y2 = 0 
            do l = l1, l2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,f,c,l) 
            enddo 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do l = l1, l2 
c              etemp = etemp + x1(a,b,l,d)*x2(e,f,c,l)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(l-l1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.6).and.(q1(3).eq.7).and.
     *    (q1(4).eq.2).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.4).and.(q2(3).eq.3).and.
     *    (q2(4).eq.7)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do l = l1, l2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,f,l,b) 
            enddo  
            enddo  
            enddo  
            
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            y2 = 0 
            do l = l1, l2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,d,c,l) 
            enddo 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do l = l1, l2 
c              etemp = etemp + x1(a,f,l,b)*x2(e,d,c,l)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(l-l1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.6).and.(q1(3).eq.3).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,f,c,k) 
            enddo 
            enddo 

         do b = b1, b2
            y1 = 0 
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,b,k,d) 
            enddo 
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,f,c,k)*x2(e,b,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.4).and.(q1(3).eq.3).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.6).and.(q2(3).eq.7).and.
     *    (q2(4).eq.2)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,d,c,k) 
            enddo 
            enddo 

         do b = b1, b2
            y1 = 0 
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,f,k,b) 
            enddo 
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,d,c,k)*x2(e,f,k,b)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.5).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.6).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
            y1 = 0 
            do b = b1, b2
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,b,e,k) 
            enddo 
            enddo 
            enddo 
         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,f,k,d) 
            enddo 
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,b,e,k)*x2(c,f,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.3).and.
     *    (q1(4).eq.7).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.6).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
         do k = k1, k2 
            y1 = y1 + 1 
            scr1(y1) = x1(a,b,c,k) 
         enddo 
         enddo 
         enddo 
         enddo 

         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,f,k,d) 
            enddo 
            y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,b,c,k)*x2(e,f,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.6).and.(q1(3).eq.7).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.2).and.(q2(3).eq.3).and.
     *    (q2(4).eq.7)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do a = a1, a2 
            do l = l1, l2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,f,l,d) 
            enddo 
            enddo 
         do c = c1, c2
         do b = b1, b2
            y2 = 0 
            do l = l1, l2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,b,c,l) 
            enddo 
            y1 = 0 
         do a = a1, a2

            etemp = 0.0 

            do l = l1, l2 
c              etemp = etemp + x1(a,f,l,d)*x2(e,b,c,l)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(l-l1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.3).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(c,k,a,d) 
            enddo 
            enddo 
      
         do b = b1, b2
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,b,k,f) 
            enddo 
            y1 = 0 
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(c,k,a,d)*x2(e,b,k,f)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.3).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.2).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.6).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
         do k = k1, k2 
            y1 = y1 + 1 
            scr1(y1) = x1(c,k,a,b) 
         enddo 
         enddo 
         enddo 
         enddo 

         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,f,k,d) 
            enddo 
            y1 = 0 
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(c,k,a,b)*x2(e,f,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.1).and.(q1(2).eq.7).and.(q1(3).eq.5).and.
     *    (q1(4).eq.4).and.  
     *    (q2(1).eq.3).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.6)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(a,k,e,d) 
            enddo 
            enddo 

         do c = c1, c2
         do b = b1, b2
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(c,b,k,f) 
            enddo 
            
            y1 = 0 
         do a = a1, a2

            etemp = 0.0 

            do k = k1, k2 
c              etemp = etemp + x1(a,k,e,d)*x2(c,b,k,f)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

      if ((q1(1).eq.3).and.(q1(2).eq.7).and.(q1(3).eq.1).and.
     *    (q1(4).eq.6).and.  
     *    (q2(1).eq.5).and.(q2(2).eq.2).and.(q2(3).eq.7).and.
     *    (q2(4).eq.4)) then   
c        print *,'cwork6441: Optimized version lineno ',current_line

         yy = 0 
         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
            y1 = 0 
            do a = a1, a2
            do k = k1, k2 
               y1 = y1 + 1 
               scr1(y1) = x1(c,k,a,f) 
            enddo 
            enddo 
         do b = b1, b2
            y2 = 0 
            do k = k1, k2 
               y2 = y2 + 1 
               scr2(y2) = x2(e,b,k,d) 
            enddo 
            y1 = 0 
         do a = a1, a2

            etemp = 0.0 ! VFL  

            do k = k1, k2 
c              etemp = etemp + x1(c,k,a,f)*x2(e,b,k,d)
               y1 = y1 + 1 
               etemp = etemp + scr1(y1)*scr2(k-k1+1)
            enddo

            yy = yy + 1 
            y(yy) = etemp  
         
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return 
c        go to 99 

      endif 

99    continue 
c-----------------------------------------------------------------------
c   Perform the contraction again in unoptimized mode.
c-----------------------------------------------------------------------

      write(6,*) 'Using unoptimized cwork6441 at line number ', 
     *            current_line 
            call cwork6441_unopt(y, na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     &                              ne1,ne2,nf1,nf2, 
     &                               a1, a2, b1, b2, c1, c2, d1, d2,
     &                               e1,e2,f1,f2, 
     &                           inda, indb, indc, indd, inde, indf,
     &                           x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     &                               m1, m2, n1, n2, r1, r2, s1, s2,
     &                           indx1,
     &                           x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     &                               i1, i2, j1, j2, k1, k2, l1, l2,
     &                           indx2,
     &                           cind, flopcount, scr1, scr2)

      return
      end


