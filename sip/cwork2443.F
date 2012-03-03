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
      subroutine cwork2443(y,na1,na2,nb1,nb2,a1,a2,b1,b2,inda, indb, 
     *                     x1,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *                     c1,c2,d1,d2,e1,e2,f1,f2, indx1,
     *                     x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                     i1,i2,j1,j2,k1,k2,l1,l2,indx2,
     *                     flopcount,scr1, scr2) 
c-------------------------------------------------------------------------
c   Performs a "2443" contraction: 
c      2 index output array
c      4 index operand array
c      4 index operand array
c      3 index contraction.
c
c   I. e., contract the 3 indices of the 2nd operand array out of the 
c   first operand array into the 2-index output array.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'
      include 'parallel_info.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,i1,i2,j1,j2,
     *        k1,k2,l1,l2
      integer inda, indb, indx1(4), indx2(4)
      integer cind1(3), cind2(3), cpos1(3), cpos2(3)
      integer flopcount

#ifdef __crayx1
      integer*8 lda, ldb, ldc, m, n, k
#else
      integer lda, ldb, ldc, m, n, k
#endif

      double precision y(na1:na2,nb1:nb2)
      double precision x1(nc1:nc2,nd1:nd2,ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2)
      double precision scr1(*), scr2(*)

      integer ia1, ia2, ib1, ib2
      integer m1,m2,m3
      integer c1start(3), c1end(3),
     *        c1dim(3), c2dim(3)
      integer ix1(4), ix2(4)
      integer ncx1
      integer i, j, l, a, b, c, d, e, f
      integer next

      double precision xval
      logical debug, flag

      debug = cwork_debug 
      flag = .false.

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y and x2 
c   indices.
c---------------------------------------------------------------------------

      ia1 = 0
      ia2 = 0
      ib1 = 0
      ib2 = 0
      do i = 1, 4
         if (indx1(i) .eq. inda) ia1 = i
         if (indx1(i) .eq. indb) ib1 = i
         if (indx2(i) .eq. inda) ia2 = i
         if (indx2(i) .eq. indb) ib2 = i
      enddo

      if (ia1 .eq. 0 .and. ib1 .eq. 0) then
         print *,'Error in cwork2443: First operand array must ',
     *           ' contain either index ',inda,' or index ',indb,'.'
         print *,'indx1 array: ',(indx1(i),i=1,4)
         call abort_job()
      endif

      if (ia2 .eq. 0 .and. ib2 .eq. 0) then
         print *,'Error in cwork2443: Second operand array must ',
     *           ' contain either index ',inda,' or index ',indb,'.'
         print *,'indx2 array: ',(indx2(i),i=1,4)
         call abort_job()
      endif

      if (ia1 .ne. 0 .and. ib1 .ne. 0) then
         print *,'Error in cwork2443: First operand array ',
     *        'contains both indices of the output.'
         print *,'indx1 array: ',(indx1(i),i=1,4)
         print *,'inda, indb = ',inda,indb
         call abort_job()
      endif

      if (ia2 .ne. 0 .and. ib2 .ne. 0) then
         print *,'Error in cwork2443: Second operand array ',
     *        'contains both indices of the output.'
         print *,'indx2 array: ',(indx2(i),i=1,4)
         print *,'inda, indb = ',inda,indb
         call abort_job()
      endif

c      if (trace) then
c         print *,'CWORK2443 at line ',current_line,' ia1,ib1 ',
c     *     ia1,ib1,' ia2,ib2, ',ia2,ib2,
c     *     ' indx1 ',(indx1(i),i=1,4),' indx2 ',(indx2(i),i=1,4)
c      endif

      if (ia1 .eq. 1 .and. ib2 .eq. 1) then
         if (indx1(2) .eq. indx2(2) .and.
     *       indx1(3) .eq. indx2(3) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(a,i,j,k) * X2(b,i,j,k)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (i2-i1+1)
            k = (d2-d1+1)*(e2-e1+1)*(f2-f1+1)
            lda = nc2-nc1+1
            ldb = (ni2-ni1+1)
            ldc = na2-na1+1
            call dgemm('N', 'T', m,n,k,1.d0, x1(c1,d1,e1,f1), lda,
     *                 x2(i1,j1,k1,l1), ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 1 .and. ib2 .eq. 4) then
         if (indx1(2) .eq. indx2(1) .and.
     *       indx1(3) .eq. indx2(2) .and.
     *       indx1(4) .eq. indx2(3)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(a,i,j,k) * X2(i,j,k,b)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (l2-l1+1)
            k = (d2-d1+1)*(e2-e1+1)*(f2-f1+1)
            lda = nc2-nc1+1
            ldb = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
            ldc = na2-na1+1
            call dgemm('N', 'N', m,n,k,1.d0, x1(c1,d1,e1,f1), lda,
     *                 x2(i1,j1,k1,l1), ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia2 .eq. 1 .and. ib1 .eq. 1) then
         if (indx1(2) .eq. indx2(2) .and.
     *       indx1(3) .eq. indx2(3) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) =  X1(b,i,j,k) * X2(a,i,j,k)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (d2-d1+1)*(e2-e1+1)*(f2-f1+1)
            lda = ni2-ni1+1
            ldb = nc2-nc1+1
            ldc = na2-na1+1
            call dgemm('N', 'T', m,n,k,1.d0, x2(i1,j1,k1,l1), lda,
     *                 x1(c1,d1,e1,f1), ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ib1 .eq. 4 .and. ia2 .eq. 4) then
         if (indx1(2) .eq. indx2(2) .and.
     *       indx1(3) .eq. indx2(3) .and.
     *       indx1(1) .eq. indx2(1)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,j,k,b) * X2(i,j,k,a)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)*(d2-d1+1)*(e2-e1+1)
            ldb = (nc2-nc1+1)*(nd2-nd1+1)*(ne2-ne1+1)
            lda = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
            ldc = na2-na1+1
            call dgemm('T', 'N', m,n,k,1.d0, x2(i1,j1,k1,l1), lda,
     *                 x1(c1,d1,e1,f1), ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ib1 .eq. 2 .and. ia2 .eq. 2) then
         if (indx1(1) .eq. indx2(1) .and.
     *       indx1(3) .eq. indx2(3) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,b,j,k) * X2(i,a,j,k)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)
            lda = (ni2-ni1+1)
            ldb = (nc2-nc1+1)
            ldc = na2-na1+1

            xval = 0.d0
            do f = f1, f2
            do e = e1, e2
               call dgemm('T', 'N', m,n,k,1.d0, x2(i1,j1,e,f), lda,
     *                 x1(c1,d1,e,f), ldb,
     *                 xval, y(a1,b1), ldc)
               xval = 1.d0
            enddo
            enddo

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 2 .and. ib2 .eq. 2) then
         if (indx1(1) .eq. indx2(1) .and.
     *       indx1(3) .eq. indx2(3) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,a,j,k) * X2(i,b,j,k)
c---------------------------------------------------------------------------
           
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)
            lda = (nc2-nc1+1)
            ldb = (ni2-ni1+1)
            ldc = na2-na1+1

            xval = 0.d0
            do f = f1, f2
            do e = e1, e2
               call dgemm('T', 'N', m,n,k,1.d0, x1(c1,d1,e,f), lda,
     *                 x2(i1,j1,e,f), ldb,
     *                 xval, y(a1,b1), ldc)
               xval = 1.d0
            enddo
            enddo

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 3 .and. ib2 .eq. 3) then
         if (indx1(1) .eq. indx2(1) .and.
     *       indx1(2) .eq. indx2(2) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,j,a,k) * X2(i,j,b,k)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)*(d2-d1+1)
            lda = (nc2-nc1+1)*(nd2-nd1+1)
            ldb = (ni2-ni1+1)*(nj2-nj1+1)
            ldc = na2-na1+1

            xval = 0.d0
            do f = f1, f2
               call dgemm('T', 'N', m,n,k,1.d0, x1(c1,d1,e1,f), lda,
     *                 x2(i1,j1,k1,f), ldb,
     *                 xval, y(a1,b1), ldc)
               xval = 1.d0
            enddo

            flag = .true.
            go to 2000
         endif
      endif

      if (ib1 .eq. 3 .and. ia2 .eq. 3) then
         if (indx1(1) .eq. indx2(1) .and.
     *       indx1(2) .eq. indx2(2) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,j,b,k) * X2(i,j,a,k)
c---------------------------------------------------------------------------

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)*(d2-d1+1)
            lda = (ni2-ni1+1)*(nj2-nj1+1)
            ldb = (nc2-nc1+1)*(nd2-nd1+1)
            ldc = na2-na1+1

            xval = 0.d0
            do f = f1, f2
               call dgemm('T', 'N', m,n,k,1.d0, x2(i1,j1,k1,f), lda,
     *                 x1(c1,d1,e1,f), ldb,
     *                 xval, y(a1,b1), ldc)
               xval = 1.d0
            enddo

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 2 .and. ib2 .eq. 2) then
         if (indx1(1) .eq. indx2(3) .and.
     *       indx1(3) .eq. indx2(1) .and.
     *       indx1(4) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,a,j,k) * X2(j,b,i,k)
c---------------------------------------------------------------------------

            next = 0
            do f = f1,f2
            do e = e1, e2
            do c = c1, c2
            do d = d1, d2
               next = next + 1
               scr1(next) = x1(c,d,e,f)
            enddo
            enddo
            enddo
            enddo
            
            call cco__remap_1234_to_2314(ni2-ni1+1,nj2-nj1+1,
     *                                   nk2-nk1+1,nl2-nl1+1,
     *                                   x2, scr2)

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)*(e2-e1+1)*(f2-f1+1)
            lda = (a2-a1+1)
            ldb = (nj2-nj1+1)
            ldc = a2-a1+1
            call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *                 scr2, ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 2 .and. ib2 .eq. 4) then
         if (indx1(1) .eq. indx2(3) .and.
     *       indx1(3) .eq. indx2(1) .and.
     *       indx1(4) .eq. indx2(2)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,a,j,k) * X2(j,k,i,b)
c---------------------------------------------------------------------------

            next = 0
            do f = f1,f2
            do e = e1, e2
            do c = c1, c2
            do d = d1, d2
               next = next + 1
               scr1(next) = x1(c,d,e,f)
            enddo
            enddo
            enddo
            enddo
            
            call cco__remap_1234_to_3124(ni2-ni1+1,nj2-nj1+1,
     *                                   nk2-nk1+1,nl2-nl1+1,
     *                                   x2, scr2)

            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)*(e2-e1+1)*(f2-f1+1)
            lda = (a2-a1+1)
            ldb = (nj2-nj1+1)*(ni2-ni1+1)*(nk2-nk1+1)
            ldc = a2-a1+1
            call dgemm('N', 'N', m,n,k,1.d0, scr1, lda,
     *                 scr2, ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 4 .and. ib2 .eq. 2) then
         if (indx1(1) .eq. indx2(1) .and.
     *       indx1(2) .eq. indx2(4) .and.
     *       indx1(3) .eq. indx2(3)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,j,k,a) * X2(i,b,k,j)
c---------------------------------------------------------------------------

            call cco__remap_1234_to_2143(ni2-ni1+1,nj2-nj1+1,
     *                                   nk2-nk1+1,nl2-nl1+1,
     *                                   x2, scr2)

           
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (c2-c1+1)*(e2-e1+1)*(d2-d1+1)
            lda = (nc2-nc1+1)*(nd2-nd1+1)*(ne2-ne1+1)
            ldb = j2-j1+1
            ldc = a2-a1+1
            call dgemm('T', 'T', m,n,k,1.d0, x1(c1,d1,e1,f1), lda,
     *                 scr2, ldb,
     *                 0.d0, y(a1,b1), ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ib1 .eq. 3 .and. ia2 .eq. 1) then
         if (indx1(1) .eq. indx2(3) .and.
     *       indx1(2) .eq. indx2(4) .and.
     *       indx1(4) .eq. indx2(2)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,j,b,k) * X2(a,k,i,j)
c---------------------------------------------------------------------------

c            call cco__remap_1234_to_1243(nc2-nc1+1,nd2-nd1+1,
c     *                                   ne2-ne1+1,nf2-nf1+1,
c     *                                   x1, scr1)

            next = 0
            do e = e1,e2
            do f = f1,f2
            do d = d1, d2
            do c = c1, c2
               next = next + 1
               scr1(next) = x1(c,d,e,f)
            enddo
            enddo
            enddo
            enddo

            next = 0
            do j = j1, j2
            do l = l1, l2
            do k = k1, k2
            do i = i1, i2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (j2-j1+1)*(k2-k1+1)*(l2-l1+1)
            lda = (i2-i1+1)
            ldb = (c2-c1+1)*(d2-d1+1)*(f2-f1+1)
            ldc = a2-a1+1

               call dgemm('N', 'N', m,n,k,1.d0, scr2, lda,
     *                 scr1, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 1 .and. ib2 .eq. 1) then
         if (indx1(2) .eq. indx2(4) .and.
     *       indx1(3) .eq. indx2(3) .and.
     *       indx1(4) .eq. indx2(2)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(a,i,j,k) * X2(b,k,j,i)
c---------------------------------------------------------------------------

            next = 0
            do j = j1, j2
            do k = k1, k2
            do l = l1, l2
            do i = i1, i2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (j2-j1+1)*(k2-k1+1)*(l2-l1+1)
            lda = (nc2-nc1+1)
            ldb = (i2-i1+1)
            ldc = a2-a1+1

               call dgemm('N', 'T', m,n,k,1.d0, x1, lda,
     *                 scr2, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 1 .and. ib2 .eq. 2) then
         if (indx1(2) .eq. indx2(1) .and.
     *       indx1(3) .eq. indx2(4) .and.
     *       indx1(4) .eq. indx2(3)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(a,i,j,k) * X2(i,b,k,j)
c---------------------------------------------------------------------------

            next = 0
            do j = j1, j2
            do k = k1, k2
            do l = l1, l2
            do i = i1, i2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (i2-i1+1)*(k2-k1+1)*(l2-l1+1)
            lda = (nc2-nc1+1)
            ldb = (i2-i1+1)*(k2-k1+1)*(l2-l1+1)
            ldc = a2-a1+1

               call dgemm('N', 'N', m,n,k,1.d0, x1, lda,
     *                 scr2, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 1 .and. ib2 .eq. 3) then
         if (indx1(2) .eq. indx2(4) .and.
     *       indx1(3) .eq. indx2(1) .and.
     *       indx1(4) .eq. indx2(2)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(a,i,j,k) * X2(j,k,b,i)
c---------------------------------------------------------------------------

            next = 0
            do k = k1, k2
            do j = j1, j2
            do i = i1, i2
            do l = l1, l2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (i2-i1+1)*(j2-j1+1)*(l2-l1+1)
            lda = (nc2-nc1+1)
            ldb = (i2-i1+1)*(j2-j1+1)*(l2-l1+1)
            ldc = a2-a1+1

               call dgemm('N', 'N', m,n,k,1.d0, x1, lda,
     *                 scr2, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 1 .and. ib2 .eq. 4) then
         if (indx1(2) .eq. indx2(3) .and.
     *       indx1(3) .eq. indx2(2) .and.
     *       indx1(4) .eq. indx2(1)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(a,i,j,k) * X2(k,j,i,b)
c---------------------------------------------------------------------------

            next = 0
            do l = l1, l2
            do i = i1, i2
            do j = j1, j2
            do k = k1, k2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
            lda = (nc2-nc1+1)
            ldb = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
            ldc = a2-a1+1

               call dgemm('N', 'N', m,n,k,1.d0, x1, lda,
     *                 scr2, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 4 .and. ib2 .eq. 3) then
         if (indx1(1) .eq. indx2(2) .and.
     *       indx1(2) .eq. indx2(1) .and.
     *       indx1(3) .eq. indx2(4)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,j,k,a) * X2(j,i,b,k)
c---------------------------------------------------------------------------

            next = 0
            do l = l1, l2
            do i = i1, i2
            do j = j1, j2
            do k = k1, k2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (i2-i1+1)*(j2-j1+1)*(l2-l1+1)
            lda = (nc2-nc1+1)*(nd2-nd1+1)*(ne2-ne1+1)
            ldb = (k2-k1+1)
            ldc = a2-a1+1

               call dgemm('T', 'T', m,n,k,1.d0, x1, lda,
     *                 scr2, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ia1 .eq. 2 .and. ib2 .eq. 1) then
         if (indx1(1) .eq. indx2(2) .and.
     *       indx1(3) .eq. indx2(4) .and.
     *       indx1(4) .eq. indx2(3)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,a,j,k) * X2(b,i,k,j)
c---------------------------------------------------------------------------

            next = 0
            do d = d1, d2
            do f = f1, f2
            do e = e1, e2
            do c = c1, c2
               next = next + 1
               scr1(next) = x1(c,d,e,f)
            enddo
            enddo
            enddo
            enddo

            next = 0
            do k = k1, k2
            do l = l1, l2
            do j = j1, j2
            do i = i1, i2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (l2-l1+1)*(j2-j1+1)*(k2-k1+1)
            lda = (c2-c1+1)*(f2-f1+1)*(e2-e1+1)
            ldb = (i2-i1+1)
            ldc = a2-a1+1

               call dgemm('T', 'T', m,n,k,1.d0, scr1, lda,
     *                 scr2, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ib1 .eq. 2 .and. ia2 .eq. 3) then
         if (indx1(1) .eq. indx2(4) .and.
     *       indx1(3) .eq. indx2(2) .and.
     *       indx1(4) .eq. indx2(1)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,b,j,k) * X2(k,j,a,i)
c---------------------------------------------------------------------------

            next = 0
            do d = d1, d2
            do f = f1, f2
            do e = e1, e2
            do c = c1, c2
               next = next + 1
               scr1(next) = x1(c,d,e,f)
            enddo
            enddo
            enddo
            enddo

            next = 0
            do i = i1, i2
            do j = j1, j2
            do l = l1, l2
            do k = k1, k2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (l2-l1+1)*(j2-j1+1)*(i2-i1+1)
            lda = (k2-k1+1)
            ldb = (c2-c1+1)*(e2-e1+1)*(f2-f1+1)
            ldc = a2-a1+1

               call dgemm('N', 'N', m,n,k,1.d0, scr2, lda,
     *                 scr1, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

      if (ib1 .eq. 2 .and. ia2 .eq. 3) then
         if (indx1(1) .eq. indx2(2) .and.
     *       indx1(3) .eq. indx2(4) .and.
     *       indx1(4) .eq. indx2(1)) then

c---------------------------------------------------------------------------
c   Y(a,b) = X1(i,b,j,k) * X2(k,i,a,j)
c---------------------------------------------------------------------------

            next = 0
            do d = d1, d2
            do f = f1, f2
            do e = e1, e2
            do c = c1, c2
               next = next + 1
               scr1(next) = x1(c,d,e,f)
            enddo
            enddo
            enddo
            enddo

            next = 0
            do i = i1, i2
            do l = l1, l2
            do j = j1, j2
            do k = k1, k2
               next = next + 1
               scr2(next) = x2(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            
            m = (a2-a1+1)
            n = (b2-b1+1)
            k = (l2-l1+1)*(j2-j1+1)*(i2-i1+1)
            lda = (k2-k1+1)
            ldb = (c2-c1+1)*(e2-e1+1)*(f2-f1+1)
            ldc = a2-a1+1

               call dgemm('N', 'N', m,n,k,1.d0, scr2, lda,
     *                 scr1, ldb,
     *                 0.d0, y, ldc)

            flag = .true.
            go to 2000
         endif
      endif

 2000 continue
      if (flag) then
         if (.not. debug) return

c--------------------------------------------------------------------------
c  Self-checking mode: Save data that was formed by the optimized contraction.
c--------------------------------------------------------------------------

         n = (b2-b1+1)*(a2-a1+1)
         call xdcopy(n, y, 1, scr1, 1)
      endif

22222 continue

c----------------------------------------------------------------------------
c  Find the contracted indices.
c----------------------------------------------------------------------------

      ncx1 = 0
      do i = 1, 4
         if (indx1(i) .ne. inda .and. 
     *       indx1(i) .ne. indb) then
           ncx1 = ncx1 + 1
           cind1(ncx1) = indx1(i)
           cpos1(ncx1) = i
           if (i .eq. 1) then
              c1start(ncx1) = c1
              c1end(ncx1)   = c2
              c1dim(ncx1)   = nc2-nc1+1
           else if (i .eq. 2) then
              c1start(ncx1) = d1
              c1end(ncx1)   = d2
              c1dim(ncx1)   = nd2-nd1+1
           else if (i .eq. 3) then
              c1start(ncx1) = e1
              c1end(ncx1)   = e2
              c1dim(ncx1)   = ne2-ne1+1
           else 
              c1start(ncx1) = f1
              c1end(ncx1)   = f2
              c1dim(ncx1)   = nf2-nf1+1
           endif 

c---------------------------------------------------------------------------
c   Find the index in indx2 corresponding to cind1(nc1).
c---------------------------------------------------------------------------

            cpos2(ncx1) = 0
            do j = 1, 4
               if (indx2(j) .eq. cind1(ncx1)) then
                  cpos2(ncx1) = j
                  cind2(ncx1) = indx2(j)
                  if (j .eq. 1) then
                     c2dim(ncx1) = i2-i1+1
                  else if (j .eq. 2) then
                     c2dim(ncx1) = j2-j1+1
                  else if (j .eq. 3) then
                     c2dim(ncx1) = k2-k1+1
                  else
                     c2dim(ncx1) = l2-l1+1
                  endif
               endif
            enddo

            if (cpos2(ncx1) .eq. 0) then
               print *,'Error in cwork2443: A contraction index ',
     *           'in the 1st operand array',
     *           ' was not found in the 2nd operand array.'
               print *,'indx2 = ',(indx2(j),j=1,4)
               print *,'indx1 = ',(indx1(j),j=1,4)
               print *,'Contraction index = ',cind1(ncx1)
               call abort_job()
            endif
         endif
      enddo

      if (ncx1 .ne. 3) then
         print *,'Error in cwork2443: First operand array must ',
     *           'contain 3 contracted indices.'
         print *,'Detected ',ncx1,' contracted indices, indices ',
     *          (cind1(i),i=1,ncx1)
         print *,'ia1,ib1,ia2,ib2 = ',ia1,ib1,ia2,ib2
         print *,'indx1: ',(indx1(i),i=1,4)
         print *,'indx2: ',(indx2(i),i=1,4)

         call abort_job()
      endif

c     write(6,*) ' Unoptimized contraction in cwork2443 line:', 
c    *             current_line 

      do b = b1,b2
         if (ib1 .ne. 0) ix1(ib1) = b 
         if (ib2 .ne. 0) ix2(ib2) = b 
      do a = a1,a2
         if (ia1 .ne. 0) ix1(ia1) = a 
         if (ia2 .ne. 0) ix2(ia2) = a 
         y(a,b) = 0.

         do m1 = c1start(1), c1end(1)
            ix1(cpos1(1)) = m1
            ix2(cpos2(1)) = m1
         do m2 = c1start(2), c1end(2)
            ix1(cpos1(2)) = m2
            ix2(cpos2(2)) = m2
         do m3 = c1start(3), c1end(3)
            ix1(cpos1(3)) = m3
            ix2(cpos2(3)) = m3
            y(a,b) = y(a,b) + x1(ix1(1),ix1(2),ix1(3),ix1(4))*
     *               x2(ix2(1),ix2(2),ix2(3),ix2(4))
         enddo
         enddo
         enddo
      enddo 
      enddo 

      if (debug .and. flag) then

c-------------------------------------------------------------------------
c   Compare the optimized results with the unoptimized result.
c-------------------------------------------------------------------------

         next = 0
         do b = b1,b2
         do a = a1,a2
            next = next + 1
            if (dabs(y(a,b) - scr1(next)) .gt. 1.d-12) then
                print *,'Error in cwork2443 contraction at line ',
     *                     current_line
                print *,'OPT: ',scr1(next),' UNOPT: ',y(a,b)
                call abort_job()
            endif
         enddo
         enddo
      endif
      return
      end
