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
      subroutine cwork2422(y,na1,na2,nb1,nb2,a1,a2,b1,b2,inda, indb, 
     *                     x1,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *                     c1,c2,d1,d2,e1,e2,f1,f2, indx1,
     *                     x2,ni1,ni2,nj1,nj2,
     *                     i1,i2,j1,j2,indi, indj,
     *                     scr1, scr2,flopcount) 
c-------------------------------------------------------------------------
c   Performs a "2422" contraction: 
c      2 index output array
c      4 index operand array
c      2 index operand array
c      2 index contraction.
c
c   I. e., contract the 2 indices of the 2nd operand array out of the 
c   first operand array into the 2-index output array.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      character*1 char(4)  
      character*2 char2a, char2b   
      character*4 char4 

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,i1,i2,j1,j2
      integer inda, indb, indx1(4), indi, indj, h(4), perm(4) 
      integer flopcount

      double precision y(na1:na2,nb1:nb2)
      double precision x1(nc1:nc2,nd1:nd2,ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2)
      double precision scr1(*), scr2(*) 

      integer ia, ib, ii, ij, ix(4)
      integer i, j, a, b, c, d, e, f 
#ifdef __crayx1
      integer*8 ncount, m, n, k, lda, ldb, ldc 
#else
      integer ncount, m, n, k, lda, ldb, ldc 
#endif

      logical check
      check = .false.

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y and x2 
c   indices.
c---------------------------------------------------------------------------

      do i = 1, 4
         if (indx1(i) .eq. inda) then 
            ia = i 
         else if (indx1(i) .eq. indb) then
            ib = i 
         else if (indx1(i) .eq. indi) then
            ii = i
         else if (indx1(i) .eq. indj) then
            ij = i
         else
            print *,'Error: Invalid index in cwork2422'
            print *,'X1 index is ',indx1(i),' y indices: ',inda, indb,
     *              ' x2 indices: ',indi,indj
            call abort_job() 
         endif
      enddo

      flopcount = 2 * (i2-i1+1)*(j2-j1+1)*(b2-b1+1)*(a2-a1+1)

c---------------------------------------------------------------------------
c   Determine permutation required to put array 2 into x2(a,b,i,j) form.  
c---------------------------------------------------------------------------

      char(1) = '1' 
      char(2) = '2' 
      char(3) = '3' 
      char(4) = '4' 

      h(1) = inda 
      h(2) = indb 
      h(3) = indi  
      h(4) = indj  

      do i = 1, 4 
      do j = 1, 4 
         if (h(i) .eq. indx1(j)) perm(i) = j 
      enddo 
      enddo 

      char4 = char(perm(1))//char(perm(2))//char(perm(3))//char(perm(4)) 

c---------------------------------------------------------------------------
c   y(a,b) = x1(a,b,i,j)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '1234') then  

         m   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(f2-f1+1)
         lda = (nc2-nc1+1)*(nd2-nd1+1)  
         ldb = k !(ni2-i1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'N', m,n,k,1.d0, x1(c1,d1,e1,f1), lda,
     *              x2(i1,j1), ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '1234' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(b,a,i,j)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '2134') then  

         m   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(f2-f1+1)
         lda = (nc2-nc1+1)*(nd2-nd1+1)  
         ldb = k ! (ni2-i1+1)*(nj2-nj1)  
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'N', m,n,k,1.d0, x1(c1,d1,e1,f1), lda,
     *              x2(i1,j1), ldb,
     *              0.d0, scr1, ldc)

         ncount = 0 
         do a = a1, a2 
         do b = b1, b2 
            ncount = ncount + 1  
            y(a,b) = scr1(ncount) 
         enddo 
         enddo 

         go to 33 
      endif ! char4 .eq. '2134' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(b,a,j,i)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '2143') then  

         ncount = 0
         do e = e1, e2
         do f = f1, f2
         do d = d1, d2
         do c = c1, c2
            ncount = ncount + 1
            scr2(ncount) = x1(c,d,e,f)
         enddo  
         enddo  
         enddo  
         enddo  

         m   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(f2-f1+1)
         lda = (c2-c1+1)*(d2-d1+1)  
         ldb = k ! (ni2-i1+1)*(nj2-nj1)  
         ldc = (a2-a1+1)*(b2-b1+1)   

         call dgemm('N', 'N', m,n,k,1.d0, scr2, lda,
     *              x2(i1,j1), ldb,
     *              0.d0, scr1, ldc)

         ncount = 0 
         do a = a1, a2 
         do b = b1, b2 
            ncount = ncount + 1  
            y(a,b) = scr1(ncount) 
         enddo 
         enddo 

         go to 33 
      endif ! char4 .eq. '2143' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(i,j,a,b)*x2(i,j).    
c        ->  x2(i,j)*x1(i,j,a,b) 
c---------------------------------------------------------------------------

      if (char4 .eq. '3412') then  

         k   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         m   = (e2-e1+1)*(f2-f1+1)
         lda = (nc2-nc1+1)*(nd2-nd1+1)  
         ldb = 1 ! (ni2-i1+1) 
         ldc = 1 ! (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'N', n,m,k,1.d0, x2(i1,j1), ldb,
     *              x1(c1,d1,e1,f1), lda,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '3412' 

      if (char4 .eq. '3421') then
c----------------------------------------------------------------------------
c   y(a,b) <-- x1(j,i,a,b) * x2(i,j)
c----------------------------------------------------------------------------

         ncount = 0
         do i = i1, i2
         do j = j1, j2
            ncount = ncount + 1
            scr2(ncount) = x2(i,j)
         enddo
         enddo
 
         k   = (c2-c1+1)*(d2-d1+1)
         m   = 1 !(i2-i1+1)
         n   = (e2-e1+1)*(f2-f1+1)
         ldb = (nc2-nc1+1)*(nd2-nd1+1)
         lda = 1 ! (ni2-i1+1)
         ldc = 1 ! (na2-na1+1)*(nb2-nb1+1)

         call dgemm('N', 'N', m,n,k,1.d0, scr2, lda,
     *              x1(c1,d1,e1,f1), ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33
      endif ! char4 .eq. '3421'

c---------------------------------------------------------------------------
c   y(a,b) = x1(a,b,j,i)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '1243') then  

         ncount = 0 
         do i = i1, i2 
         do j = j1, j2 
            ncount = ncount + 1 
            scr1(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)*(j2-j1+1)
         k   = (e2-e1+1)*(f2-f1+1)
         lda = (nc2-nc1+1)*(nd2-nd1+1)  
         ldb = 1 !(ni2-i1+1)*(nj2-nj1+1)  
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, x1(c1,d1,e1,f1), lda,
     *              scr1, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '1243' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(b,i,j,a)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '4123') then  

         ncount = 0
         do e = e1, e2
         do d = d1, d2
         do f = f1, f2
         do c = c1, c2
            ncount = ncount + 1
            scr1(ncount) = x1(c,d,e,f) 
         enddo
         enddo
         enddo
         enddo
 
         k   = (d2-d1+1)*(e2-e1+1)
         n   = 1 
         m   = (c2-c1+1)*(f2-f1+1)
         lda = (c2-c1+1)*(f2-f1+1)
         ldb = (ni2-ni1+1)*(nj2-nj1+1)
         ldc = (a2-a1+1)*(b2-b1+1)

         call dgemm('N', 'N', m,n,k,1.d0, scr1, lda,
     *              x2(i1,j1), ldb,
     *              0.d0, scr2, ldc)

         ncount = 0
         do a = a1, a2
         do b = b1, b2
            ncount = ncount + 1
            y(a,b) = scr2(ncount)
         enddo
         enddo

         go to 33 
      endif ! char4 .eq. '4123' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(i,j,b,a)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '4312') then  

         ncount = 0
         do e = e1, e2
         do f = f1, f2
         do d = d1, d2
         do c = c1, c2
            ncount = ncount + 1
            scr1(ncount) = x1(c,d,e,f) 
         enddo
         enddo
         enddo
         enddo
 
         k   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         m   = (e2-e1+1)*(f2-f1+1)
         lda = (c2-c1+1)*(d2-d1+1)
         ldb = 1 ! (ni2-i1+1)
         ldc = 1 ! (na2-na1+1)*(nb2-nb1+1)

         call dgemm('N', 'N', n,m,k,1.d0, x2(i1,j1), ldb,
     *              scr1, lda,
     *              0.d0, y(a1,b1), ldc)
         go to 33 
      endif ! char4 .eq. '4312' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(j,i,b,a)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '4321') then  

         ncount = 0
         do e = e1, e2
         do f = f1, f2
         do d = d1, d2
         do c = c1, c2
            ncount = ncount + 1
            scr1(ncount) = x1(c,d,e,f) 
         enddo
         enddo
         enddo
         enddo
 
         ncount = 0
         do i = i1, i2
         do j = j1, j2
            ncount = ncount + 1
            scr2(ncount) = x2(i,j)
         enddo
         enddo

         k   = (c2-c1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         m   = (e2-e1+1)*(f2-f1+1)
         lda = (c2-c1+1)*(d2-d1+1)
         ldb = (i2-i1+1)*(j2-j1+1)
         ldc = (na2-na1+1)*(nb2-nb1+1)

         call dgemm('T', 'N', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)
         go to 33 
      endif ! char4 .eq. '4321' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(a,i,b,j)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '1324') then  

         ncount = 0 
         do f = f1, f2 
         do d = d1, d2 
         do e = e1, e2 
         do c = c1, c2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (c2-c1+1)*(e2-e1+1)
         n   = 1 !(i2-i1+1)
         k   = (d2-d1+1)*(f2-f1+1)
         lda = (nc2-nc1+1)*(ne2-ne1+1)  
         ldb = 1 !(ni2-i1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '1324' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(a,i,j,b)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '1342') then  

         ncount = 0 
         do d = d1, d2 
         do f = f1, f2 
         do e = e1, e2 
         do c = c1, c2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (c2-c1+1)*(e2-e1+1)
         n   = 1 !(i2-i1+1)
         k   = (d2-d1+1)*(f2-f1+1)
         lda = (nc2-nc1+1)*(ne2-ne1+1)  
         ldb = 1 !(ni2-i1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '1342' 


c---------------------------------------------------------------------------
c   y(a,b) = x1(a,i,j,b)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '1423') then  

         ncount = 0 
         do e = e1, e2 
         do d = d1, d2 
         do f = f1, f2 
         do c = c1, c2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (c2-c1+1)*(f2-f1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(d2-d1+1)
         lda = (nc2-nc1+1)*(nf2-nf1+1)  
         ldb = 1 ! (ni2-i1+1)*(nj2-nj1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '1423' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(a,j,i,b)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '1432') then  

         ncount = 0 
         do d = d1, d2 
         do e = e1, e2 
         do f = f1, f2 
         do c = c1, c2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (c2-c1+1)*(f2-f1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(d2-d1+1)
         lda = (nc2-nc1+1)*(nf2-nf1+1)  
         ldb = 1 ! (ni2-i1+1)*(nj2-nj1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '1432' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(i,a,b,j)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '2314') then  

         ncount = 0 
         do f = f1, f2 
         do c = c1, c2 
         do e = e1, e2 
         do d = d1, d2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (e2-e1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (f2-f1+1)*(c2-c1+1)
         lda = (ne2-ne1+1)*(nd2-nd1+1)  
         ldb = 1 !(ni2-i1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '2314' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(i,a,b,j)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '2341') then  

         ncount = 0 
         do c = c1, c2 
         do f = f1, f2 
         do e = e1, e2 
         do d = d1, d2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (e2-e1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (f2-f1+1)*(c2-c1+1)
         lda = (ne2-ne1+1)*(nd2-nd1+1)  
         ldb = 1 !(ni2-i1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '2341' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(i,a,j,b)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '2413') then  
         ncount = 0 
         do e = e1, e2 
         do c = c1, c2 
         do f = f1, f2 
         do d = d1, d2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (f2-f1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(c2-c1+1)
         lda = (nf2-nf1+1)*(nd2-nd1+1)  
         ldb = 1 ! (ni2-i1+1)*(nj2-nj1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '2413' 

c---------------------------------------------------------------------------
c   y(a,b) = x1(b,j,i,b)*x2(i,j).    
c---------------------------------------------------------------------------

      if (char4 .eq. '2431') then  
         ncount = 0 
         do c = c1, c2 
         do e = e1, e2 
         do f = f1, f2 
         do d = d1, d2 
            ncount = ncount + 1 
            scr1(ncount) = x1(c,d,e,f) 
         enddo 
         enddo 
         enddo 
         enddo 

         ncount = 0 
         do j = j1, j2 
         do i = i1, i2 
            ncount = ncount + 1 
            scr2(ncount) = x2(i,j) 
         enddo 
         enddo 

         m   = (f2-f1+1)*(d2-d1+1)
         n   = 1 !(i2-i1+1)
         k   = (e2-e1+1)*(c2-c1+1)
         lda = (nf2-nf1+1)*(nd2-nd1+1)  
         ldb = 1 ! (ni2-i1+1)*(nj2-nj1+1) 
         ldc = (na2-na1+1)*(nb2-nb1+1)   

         call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *              scr2, ldb,
     *              0.d0, y(a1,b1), ldc)

         go to 33 
      endif ! char4 .eq. '2431' 

      go to 10101

33    continue 

      if (check) then
         ncount = 0
         do b = b1, b2
         do a = a1, a2
           ncount = ncount + 1
           scr2(ncount) = y(a,b)
         enddo
         enddo
      else
         return
      endif

10101 continue

      do a = a1, a2
         ix(ia) = a
      do b = b1, b2
         ix(ib) = b
         y(a,b) = 0.
 
         do j = j1, j2
            ix(ij) = j
         do i = i1, i2 
            ix(ii) = i
            y(a,b) = y(a,b) + x2(i,j)*x1(ix(1),ix(2),ix(3),ix(4))
         enddo
         enddo
      enddo
      enddo

      if (check) then
         ncount = 0
         do b = b1, b2
         do a = a1, a2
            ncount = ncount + 1
            if (dabs(scr2(ncount) - y(a,b)) .gt. 1.d-6) then
               print *,'Error in contraction at line ',current_line
               print *,'c1,c2,d1,d2,e1,e2,f1,f2 ',
     *             c1,c2,d1,d2,e1,e2,f1,f2,' i1, i2, j1, j2 ',
     *             i1, i2, j1, j2
               print *,'char4 ',char4
               call abort_job()
            endif
         enddo
         enddo
      endif

c       write(6,*) ' PERMUTATION ', char4, 
c     *     ' NOT OPTIMIZED IN cwork2422 line ',current_line 

      return
      end

