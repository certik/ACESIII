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
       subroutine cwork2644(y, na1,na2,nb1,nb2,a1,a2,b1,b2,
     *                      inda, indb,
c
     *                      x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *                         nc1,nc2,nd1,nd2,
     *                          m1, m2, n1, n2, r1, r2, s1, s2,
     *                          c1, c2, d1, d2,
     *                      indx1,
c
     *                      x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                          i1, i2, j1, j2, k1, k2, l1, l2,
     *                      indx2,
     *                      cind, flopcount, scr1, scr2, scr3)
c-------------------------------------------------------------------------
c   Performs a "6441" contraction:
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
      integer inda, indb, indc, indd, inde, indf, indx1(6), indx2(4)
      integer cind(4), flopcount

      double precision y(na1:na2,nb1:nb2) 
c     double precision y(*)
      double precision x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,
     *                    nc1:nc2,nd1:nd2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
      double precision scr1(*), scr2(*), scr3(*)
      double precision etemp  
c
c Temporary arrays used in array permutations 
c 
c END Temporary arrays used in array permutations 

c     integer p1, p2, q1, q2
      integer na, nb, nc, nd, ne, nf, ng, nh, ni, nj, nk, nl, np, nq
      integer p0(2), q1(6), q2(4) 

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
     &    c1.ne.nc1.or.c2.ne.nc2.or.d1.ne.nd1.or.d2.ne.nd2   ) then
         print *, '@CWORK2644: Assertion failed.'
         print *, '            The x1 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert x2 is contiguous
      if (i1.ne.ni1.or.j1.ne.nj1.or.k1.ne.nk1.or.l1.ne.nl1.or.
     &    i2.ne.ni2.or.j2.ne.nj2.or.k2.ne.nk2.or.l2.ne.nl2    ) then
         print *, '@CWORK2644: Assertion failed.'
         print *, '            The x2 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert y is contiguous
      if (a1.ne.na1.or.b1.ne.nb1) then 
         print *, '@CWORK2644: Assertion failed.'
         print *, '            The product tensor is not contiguous.'
         call abort_job()
      end if

c     end if (.false.)
      end if

c     go to 99 
c Try optimized version !!! VFL 

      do i = 1, 2 
         p0(i) = 0
      enddo  
      do i = 1, 6 
         q1(i) = i
      enddo
      do i = 1, 4 
         q2(i) = 0
      enddo
      
      do i = 1, 6 
         if (indx1(i) .eq. inda) p0(1) = q1(i) 
         if (indx1(i) .eq. indb) p0(2) = q1(i) 
      enddo
      
      do i = 1, 6
      do j = 1, 4 
         if (indx1(i) .eq. indx2(j)) then
            q2(j) = q1(i)
         endif
      enddo
      enddo

      if (q1(1).eq.1.and.q1(2).eq.2.and.q1(3).eq.3.and.q1(4).eq.4.and.
     *    q1(5).eq.5.and.q1(6).eq.6.and. 
     *    q2(1).eq.3.and.q2(2).eq.4.and.q2(3).eq.5.and.q2(4).eq.6) 
     *    then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do n = n1, n2
         do m = m1, m2
            y(m,n) = 0.0  
         enddo
         enddo

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(m,n) = y(m,n) + x1(m,n,r,s,c,d)*x2(r,s,c,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.1.and.q2(2).eq.2.and.q2(3).eq.3.and.q2(4).eq.4  
     *      .and.p0(1).eq.5.and.p0(2).eq.6) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do d = d1, d2
         do c = c1, c2
            etemp = 0.0 
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            etemp = etemp + x1(m,n,r,s,c,d)*x2(m,n,r,s)

         enddo
         enddo
         enddo
         enddo
            y(c,d) = etemp 
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.1.and.q2(2).eq.2.and.q2(3).eq.3.and.q2(4).eq.4  
     *      .and.p0(1).eq.6.and.p0(2).eq.5) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do d = d1, d2
         do c = c1, c2
            etemp = 0.0 
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            etemp = etemp + x1(m,n,r,s,c,d)*x2(m,n,r,s)

         enddo
         enddo
         enddo
         enddo
            y(d,c) = etemp 
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.3.and.q2(2).eq.2.and.q2(3).eq.5.and.q2(4).eq.6  
     *      .and.p0(1).eq.1.and.p0(2).eq.4) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do m = m1, m2
         do s = s1, s2
            y(m,s) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(m,s) = y(m,s) + x1(m,n,r,s,c,d)*x2(r,n,c,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.1.and.q2(2).eq.4.and.q2(3).eq.5.and.q2(4).eq.6  
     *      .and.p0(1).eq.3.and.p0(2).eq.2) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do r = r1, r2
         do n = n1, n2
            y(r,n) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(r,n) = y(r,n) + x1(m,n,r,s,c,d)*x2(m,s,c,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.1.and.q2(2).eq.2.and.q2(3).eq.5.and.q2(4).eq.6  
     *      .and.p0(1).eq.3.and.p0(2).eq.4) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do r = r1, r2
         do s = s1, s2
            y(r,s) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(r,s) = y(r,s) + x1(m,n,r,s,c,d)*x2(m,n,c,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.3.and.q2(2).eq.4.and.q2(3).eq.1.and.q2(4).eq.6  
     *      .and.p0(1).eq.5.and.p0(2).eq.2) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do c = c1, c2
         do n = n1, n2
            y(c,n) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(c,n) = y(c,n) + x1(m,n,r,s,c,d)*x2(r,s,m,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.3.and.q2(2).eq.4.and.q2(3).eq.5.and.q2(4).eq.2  
     *      .and.p0(1).eq.1.and.p0(2).eq.6) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do m = m1, m2
         do d = d1, d2
            y(m,d) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(m,d) = y(m,d) + x1(m,n,r,s,c,d)*x2(r,s,c,n)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.1.and.q2(2).eq.2.and.q2(3).eq.5.and.q2(4).eq.6  
     *      .and.p0(1).eq.3.and.p0(2).eq.4) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do r = r1, r2
         do s = s1, s2
            y(r,s) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(r,s) = y(r,s) + x1(m,n,r,s,c,d)*x2(m,n,c,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.3.and.q2(2).eq.2.and.q2(3).eq.1.and.q2(4).eq.6  
     *      .and.p0(1).eq.5.and.p0(2).eq.4) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do c = c1, c2
         do s = s1, s2
            y(c,s) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(c,s) = y(c,s) + x1(m,n,r,s,c,d)*x2(r,n,m,d)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

      if (q2(1).eq.1.and.q2(2).eq.4.and.q2(3).eq.5.and.q2(4).eq.2  
     *      .and.p0(1).eq.3.and.p0(2).eq.6) then 
c        print *,'cwork2644: Optimized version lineno ',current_line 

         do r = r1, r2
         do d = d1, d2
            y(r,d) = 0.0 
         enddo 
         enddo 

         do d = d1, d2
         do c = c1, c2
         do s = s1, s2
         do r = r1, r2
         do n = n1, n2
         do m = m1, m2

            y(r,d) = y(r,d) + x1(m,n,r,s,c,d)*x2(m,s,c,n)

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         return

      endif   

99    continue 
c-----------------------------------------------------------------------
c   Perform the contraction again in unoptimized mode.
c-----------------------------------------------------------------------

            call cwork2644_unopt(y, na1,na2,nb1,nb2, 
     &                               a1, a2, b1, b2, 
     &                           inda, indb, 
     &                           x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     &                               nc1,nc2,nd1,nd2,
     &                               m1, m2, n1, n2, r1, r2, s1, s2,
     &                               c1, c2, d1, d2,
     &                           indx1,
     &                           x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     &                               i1, i2, j1, j2, k1, k2, l1, l2,
     &                           indx2,
     &                           cind, flopcount, scr1, scr2)

      return
      end


