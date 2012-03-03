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
      subroutine cwork4421(y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                     a1,a2,b1,b2,c1,c2,d1,d2,inda, indb, 
     *                     indc, indd, x1,ne1,ne2,nf1,nf2,ng1,ng2,
     *                     nh1,nh2,e1,e2,f1,f2,g1,g2,h1,h2, 
     *                     indx1, x2,ni1,ni2,nj1,nj2,
     *                     i1,i2,j1,j2,indi, indj,
     *                     cind, flopcount,scr1, scr2) 
c-------------------------------------------------------------------------
c   Performs a "4421" contraction: 
c      4 index output array
c      4 index operand array
c      2 index operand array
c      1 index contraction.
c
c   I. e., contract one of the the 2 indices of the 2nd operand array out 
c   of the first operand array, replacing the index by the 2nd operand's 
c   2nd index.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2
      integer ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2
      integer ni1, ni2, nj1, nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,h1,h2,
     *        i1,i2,j1,j2
      integer inda, indb, indc, indd, indx1(4), indi, indj
      integer cind(2), flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2)
      double precision x2(ni1:ni2,nj1:nj2)

      double precision scr1(*)
      double precision scr2(*)
      integer next

      integer ia, ib, ic, id, icontract, ix(4), iy
      integer i, j, a, b, c, d
      integer indy(4), map(4)
      integer ie, if, ig, ih 
      integer e, f, g, h
      integer alast, dlast
      integer iprint
      double precision cval

      integer*8 ixx, c_loc64

#ifdef __crayx1
      integer*8 m, n, k, lda, ldb, ldc
#else
      integer m, n, k, lda, ldb, ldc
#endif

      iprint = cwork_iprint

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
            icontract = i
         else
            print *,'Error: Invalid index in cwork4421'
            print *,'X1 index is ',indx1(i),' y indices: ',inda, indb,
     *              indc, indd,
     *              ' x2 indices: ',indi,indj,' cind(1) = ',cind(1)
            call abort_job() 
         endif
      enddo

      indy(1) = ia
      indy(2) = ib
      indy(3) = ic
      indy(4) = id
      do i = 1, 4
         if (indy(i) .eq. 0) indy(i) = icontract
      enddo

      if (indi .eq. cind(1)) then

c--------------------------------------------------------------------------
c   Contract over the first index in X2.
c--------------------------------------------------------------------------

         flopcount = 2*(i2-i1+1)*(d2-d1+1)*(c2-c1+1)*(b2-b1+1)*
     *                 (a2-a1+1) 
         if (indy(1) .eq. 1 .and. indy(2) .eq. 2 .and. 
     *       indy(3) .eq. 3 .and.
     *       indy(4) .eq. 4) then 
            
            if (icontract .eq. 1) then 

c--------------------------------------------------------------------------
c   Contract over the 1st index of X1.
c
c   Y(a,b,c,d) = X1(i,b,c,d)*(X2(i,a)
c--------------------------------------------------------------------------

               m = (a2-a1+1)
               n = (b2-b1+1)*(c2-c1+1)*(d2-d1+1)
               k = i2-i1+1
               lda = (ni2-ni1+1)
               ldb = ne2-ne1+1
               ldc = (na2-na1+1)
               call dgemm('T', 'N', m,n,k,1.d0, x2(i1,j1), lda,
     *                 x1(i1,b1,c1,d1), ldb,
     *                 0.d0, y(a1,b1,c1,d1), ldc)
               return

            endif

            if (icontract .eq. 2) then 

c--------------------------------------------------------------------------
c   Contract over the 2nd index of X1.
c
c   Y(a,b,c,d) = X1(a,i,c,d) * X2(i,b)
c--------------------------------------------------------------------------

               m = (a2-a1+1)
               n = (b2-b1+1)
               k = i2-i1+1
               lda = (ne2-ne1+1)
               ldb = (ni2-ni1+1)
               ldc = (na2-na1+1)
               do c = c1, c2
               do d = d1, d2
                  call dgemm('N', 'N', m,n,k,1.d0, x1(a1,i1,c,d), lda,
     *                        x2(i1,j1), ldb, 0.d0,
     *                        y(a1,b1,c,d), ldc)
               enddo
               enddo

               return
            endif

            if (icontract .eq. 3) then 

c--------------------------------------------------------------------------
c   Contract over the 3rd index of X1.
c
c   Y(a,b,c,d) + X1(a,b,i,d)*X2(i,c)
c--------------------------------------------------------------------------

               m = (a2-a1+1)*(b2-b1+1)
               n = (c2-c1+1)
               k = i2-i1+1
               lda = (ne2-ne1+1)*(nf2-nf1+1)
               ldb = (ni2-ni1+1)
               ldc = (na2-na1+1)*(nb2-nb1+1)
               do d = d1, d2
                  call dgemm('N', 'N', m,n,k,1.d0, x1(a1,b1,i1,d), lda,
     *                        x2(i1,j1), ldb, 0.d0,
     *                        y(a1,b1,c1,d), ldc)
               enddo
               return
            endif

            if (icontract .eq. 4) then 

c--------------------------------------------------------------------------
c   Contract over the 4th index of X1.
c
c   Y(a,b,c,d) + X1(a,b,c,i)*X2(i,d)
c--------------------------------------------------------------------------

               m = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               n = d2-d1+1
               k = i2-i1+1
               lda = (ne2-ne1+1)*(nf2-nf1+1)*(ng2-ng1+1)
               ldb = ni2-ni1+1
               ldc = (na2-na1+1)*(nb2-nb1+1)*(nc2-nc1+1)
               call dgemm('N', 'N', m,n,k,1.d0, x1(a1,b1,c1,i1),
     *                 lda, x2(i1,j1), ldb,
     *                 0.d0, y(a1,b1,c1,d1), ldc)
               return
            endif
         else

            if (ia .eq. 2 .and. ic .eq. 3 .and. id .eq. 4 .and. 
     *          icontract .eq. 1) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(i,a,c,d) * X2(i,b)
c----------------------------------------------------------------------------

               m = (a2-a1+1)
               n = b2-b1+1
               k = i2-i1+1
               lda = (ne2-ne1+1)
               ldb = ni2-ni1+1
               ldc = (na2-na1+1)

               do d = d1, d2
               do c = c1, c2
                  call dgemm('T', 'N', m,n,k,1.d0, x1(i1,a1,c,d),
     *                 lda, x2(i1,j1), ldb,
     *                 0.d0, y(a1,b1,c,d), ldc)
               enddo
               enddo

               return
            endif

            if (ia .eq. 1 .and. ib .eq. 2 .and. ic .eq. 4 .and. 
     *          icontract .eq. 3) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(a,b,i,c) * X2(i,d)
c----------------------------------------------------------------------------

               call cco__remap_1234_to_1243(ne2-ne1+1,nf2-nf1+1,
     *                                      ng2-ng1+1,nh2-nh1+1,
     *                                      x1, scr1)
               m = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               n = d2-d1+1
               k = i2-i1+1
               lda = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               ldb = ni2-ni1+1
               ldc = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               call dgemm('N', 'N', m,n,k,1.d0, scr1,
     *                 lda, x2(i1,d1), ldb,
     *                 0.d0, y, ldc)

               return
            endif

            if (ia .eq. 4 .and. ib .eq. 1 .and. ic .eq. 2) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(b,c,i,a) * X2(i,d)
c----------------------------------------------------------------------------

               m = (b2-b1+1)*(c2-c1+1)
               n = d2-d1+1
               k = i2-i1+1
               lda = (ne2-ne1+1)*(nf2-nf1+1)
               ldb = ni2-ni1+1
               ldc = (b2-b1+1)*(c2-c1+1)

               do a = a1, a2
                  call dgemm('N', 'N', m,n,k,1.d0, x1(e1,f1,g1,a),
     *                 lda, x2(i1,d1), ldb,
     *                 0.d0, scr2, ldc)

                  next = 0
                  do d = d1,d2
                  do c = c1,c2
                  do b = b1,b2
                     next = next + 1
                     y(a,b,c,d) = scr2(next)
                  enddo
                  enddo
                  enddo
               enddo

               return
            endif

            if (ia .eq. 1 .and. ib .eq. 4 .and. ic .eq. 3) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(a,i,c,b) * X2(i,d)
c----------------------------------------------------------------------------

               call cco__remap_1234_to_1432(ne2-ne1+1,nf2-nf1+1,
     *                                      ng2-ng1+1,nh2-nh1+1,
     *                                      x1, scr1)

               m = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               n = d2-d1+1
               k = i2-i1+1
               lda = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               ldb = ni2-ni1+1
               ldc = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)

               call dgemm('N', 'N', m,n,k,1.d0, scr1,
     *                 lda, x2(i1,d1), ldb,
     *                 0.d0, y, ldc)

               return
            endif

            if (ia .eq. 3 .and. ib .eq. 4 .and. ic .eq. 1) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(c,i,a,b) * X2(i,d)
c----------------------------------------------------------------------------

               next = 0
               do f = f1,f2
               do e = e1,e2
               do h = h1, h2
               do g = g1, g2
                  next =next + 1
                  scr1(next) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo

               m = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               n = d2-d1+1
               k = i2-i1+1
               lda = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               ldb = ni2-ni1+1
               ldc = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               call dgemm('N', 'N', m,n,k,1.d0, scr1,
     *                 lda, x2(i1,d1), ldb,
     *                 0.d0, y(a1,b1,c1,d1), ldc)

               return
            endif

            if (id .eq. 4 .and. ib .eq. 3 .and. ic .eq. 1) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(c,i,b,d) * X2(i,a)
c----------------------------------------------------------------------------

               next = 0
               do f = f1, f2
               do h = h1, h2
               do g = g1, g2
               do e = e1, e2
                  next = next + 1
                  scr1(next) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo

               m = (a2-a1+1)
               n = (c2-c1+1)*(b2-b1+1)*(d2-d1+1)
               k = i2-i1+1
               lda = (ni2-ni1+1)
               ldb = (c2-c1+1)*(b2-b1+1)*(d2-d1+1)
               ldc = (a2-a1+1)

               call dgemm('T', 'T', m,n,k,1.d0, x2(i1,a1),
     *                 lda, scr1, ldb,
     *                 0.d0, scr2, ldc)

               next = 0
               do d = d1,d2
               do b = b1, b2
               do c = c1,c2
               do a = a1,a2
                  next = next + 1
                  y(a,b,c,d) = scr2(next)
               enddo
               enddo
               enddo
               enddo

               return
            endif

            if (id .eq. 1 .and. ic .eq. 2 .and. ia .eq. 4) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(d,c,i,a) * X2(i,b)
c----------------------------------------------------------------------------

               next = 0
               do g = g1, g2
               do e = e1, e2
               do f = f1, f2
               do h = h1, h2
                  next = next + 1
                  scr1(next) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo

               m = (a2-a1+1)*(c2-c1+1)*(d2-d1+1)
               n = (b2-b1+1)
               k = i2-i1+1
               lda = (a2-a1+1)*(c2-c1+1)*(d2-d1+1)
               ldb = (ni2-ni1+1)
               ldc = (a2-a1+1)*(c2-c1+1)*(d2-d1+1)
               call dgemm('N', 'N', m,n,k,1.d0, scr1,
     *                 lda, x2(i1,b1), ldb,
     *                 0.d0, scr2, ldc)

               next = 0
               do b = b1, b2
               do d = d1,d2
               do c = c1,c2
               do a = a1,a2
                  next = next + 1
                  y(a,b,c,d) = scr2(next)
               enddo
               enddo
               enddo
               enddo

               return
            endif

c----------------------------------------------------------------------------
c   Indirect indexing.
c----------------------------------------------------------------------------

c        write(6,*) ' Entering unoptimized contraction in cwork4421.F'
c    *             ,'at line number', current_line, 'cind = ', cind(1)


            do a = a1, a2
               ix(ia) = a
               if (indj .eq. inda) iy = a
            do b = b1, b2
               ix(ib) = b
               if (indj .eq. indb) iy = b
            do c = c1, c2
               ix(ic) = c
               if (indj .eq. indc) iy = c
            do d = d1, d2
               ix(id) = d
               if (indj .eq. indd) iy = d
               y(a,b,c,d) = 0.
   
               do i = i1, i2 
                  ix(icontract) = i
                  y(a,b,c,d) = y(a,b,c,d) + 
     *                x2(i,iy)*x1(ix(1),ix(2),ix(3),ix(4))
               enddo
            enddo
            enddo
            enddo
            enddo

         endif
      else

c---------------------------------------------------------------------------
c   Contract over the 2nd index in X2.
c---------------------------------------------------------------------------

         flopcount = 2*(j2-j1+1)*(d2-d1+1)*(c2-c1+1)*(b2-b1+1)*
     *                 (a2-a1+1) 
         if (indy(1) .eq. 1 .and. indy(2) .eq. 2 .and. 
     *       indy(3) .eq. 3 .and.
     *       indy(4) .eq. 4) go to 2000

         if (ia .eq. 2 .and. ib .eq. 1 .and. id .eq. 3) then 

c--------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(b,a,d,j)*X2(c,j)
c--------------------------------------------------------------------------

            m = (b2-b1+1)*(a2-a1+1)*(d2-d1+1)
            n = (c2-c1+1)
            k = j2-j1+1
            lda = (ne2-ne1+1)*(nf2-nf1+1)*(ng2-ng1+1)
            ldb = (ni2-ni1+1)
            ldc = (b2-b1+1)*(a2-a1+1)*(d2-d1+1)
            call dgemm('N', 'T', m,n,k,1.d0, x1(e1,f1,g1,h1),
     *                 lda, x2(i1,j1), ldb,
     *                 0.d0, scr1, ldc)
 
            next = 0
            do c = c1,c2
            do d = d1, d2
            do a = a1,a2
            do b = b1, b2
               next = next + 1
               y(a,b,c,d) = scr1(next)
            enddo
            enddo
            enddo
            enddo

            return
         endif
         
         if (id .eq. 4 .and. ib .eq. 2 .and. ic .eq. 1) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(c,b,j,d) * X2(a,j)
c----------------------------------------------------------------------------

               next = 0
               do g = g1, g2
               do h = h1, h2
               do e = e1, e2
               do f = f1, f2
                  next = next + 1
                  scr1(next) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo
               
               m = (a2-a1+1)
               n = (b2-b1+1)*(c2-c1+1)*(d2-d1+1)
               k = j2-j1+1
               lda = ni2-ni1+1             
               ldb = (b2-b1+1)*(c2-c1+1)*(d2-d1+1)
               ldc = (a2-a1+1)
               call dgemm('N', 'T', m,n,k,1.d0, x2(a1,j1),
     *                 lda, scr1, ldb,
     *                 0.d0, y, ldc)

               return
            endif

         if (id .eq. 4 .and. ib .eq. 1 .and. ic .eq. 3) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(b,j,c,d) * X2(a,j)
c----------------------------------------------------------------------------

               next = 0
               do h = h1, h2
               do g = g1, g2
               do e = e1, e2
               do f = f1, f2
                  next = next + 1
                  scr1(next) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo
               
               m = (a2-a1+1)
               n = (b2-b1+1)*(c2-c1+1)*(d2-d1+1)
               k = j2-j1+1
               lda = ni2-ni1+1             
               ldb = f2-f1+1
               ldc = (a2-a1+1)
               call dgemm('N', 'N', m,n,k,1.d0, x2(a1,j1),
     *                 lda, scr1, ldb,
     *                 0.d0, y, ldc)

               return
            endif

            if (id .eq. 1 .and. ib .eq. 4 .and. ic .eq. 3) then

c----------------------------------------------------------------------------
c   Y(a,b,c,d) = X1(d,j,c,b) * X2(a,j)
c----------------------------------------------------------------------------

               call cco__remap_1234_to_4312(ne2-ne1+1,nf2-nf1+1,
     *                                     ng2-ng1+1,nh2-nh1+1,
     *                                     x1, scr1)         

               m = (a2-a1+1)
               n = (b2-b1+1)*(c2-c1+1)*(d2-d1+1)
               k = j2-j1+1
               lda = ni2-ni1+1
               ldb = (nh2-nh1+1)*(ng2-ng1+1)*(ne2-ne1+1)
               ldc = (a2-a1+1)
               call dgemm('N', 'T', m,n,k,1.d0, x2(a1,j1),
     *                 lda, scr1, ldb,
     *                 0.d0, y, ldc)

               return
            endif

c----------------------------------------------------------------------------
c   Indirect indexing.
c----------------------------------------------------------------------------

c        write(6,*) ' Entering unoptimized contraction in cwork4421.F'
c    *             ,'at line number', current_line, 'cind = ', cind(1)


         do a = a1, a2
            if (ia .eq. 0) then
               do j = j1, j2
                  scr1(j) = x2(a,j)
               enddo
            else
               ix(ia) = a
            endif

         do b = b1, b2
            if (ib .eq. 0) then
                do j = j1, j2
                  scr1(j) = x2(b,j)
               enddo
            else
               ix(ib) = b
            endif

         do c = c1, c2
            if (ic .eq. 0) then
                do j = j1, j2
                  scr1(j) = x2(c,j)
               enddo
            else
               ix(ic) = c
            endif

         do d = d1, d2
            if (id .eq. 0) then
                do j = j1, j2
                  scr1(j) = x2(d,j)
               enddo
            else
               ix(id) = d
            endif
            y(a,b,c,d) = 0.

            do j = j1, j2 
               ix(icontract) = j
               y(a,b,c,d) = y(a,b,c,d) + 
     *                scr1(j)*x1(ix(1),ix(2),ix(3),ix(4))
            enddo
         enddo
         enddo
         enddo
         enddo

         return

 2000    continue

c---------------------------------------------------------------------------
c   Optimized contractions.
c---------------------------------------------------------------------------

         if (icontract .eq. 1) then 

c--------------------------------------------------------------------------
c   Contract over the 1st index of X1.
c
c   Y(a,b,c,d) = X1(j,b,c,d) * X2(a,j)
c--------------------------------------------------------------------------

               m = (a2-a1+1)
               n = (b2-b1+1)*(c2-c1+1)*(d2-d1+1)
               k = j2-j1+1
               lda = (ni2-ni1+1)
               ldb = ne2-ne1+1
               ldc = (na2-na1+1)
               call dgemm('N', 'N', m,n,k,1.d0, x2(i1,j1),
     *                 lda, x1(j1,b1,c1,d1), ldb,
     *                 0.d0, y(a1,b1,c1,d1), ldc)
            return
         endif

         if (icontract .eq. 2) then 

c--------------------------------------------------------------------------
c   Contract over the 2nd index of X1.
c
c   Y(a,b,c,d) = X1(a,j,c,d) * X2(b,j)
c--------------------------------------------------------------------------

               m = (a2-a1+1)
               n = (b2-b1+1)
               k = j2-j1+1
               lda = (ne2-ne1+1)
               ldb = (ni2-ni1+1)
               ldc = (na2-na1+1)
               do c = c1, c2
               do d = d1, d2
                  call dgemm('N', 'T', m,n,k,1.d0, x1(a1,j1,c,d), lda,
     *                        x2(i1,j1), ldb, 0.d0,
     *                        y(a1,b1,c,d), ldc)
               enddo
               enddo
               return
         endif

         if (icontract .eq. 3) then 

c--------------------------------------------------------------------------
c   Contract over the 3rd index of X1.
c
c   Y(a,b,c,d) = X1(a,b,j,d) * X2(c,j)
c--------------------------------------------------------------------------

c--------------------------------------------------------------------------
c   First, sort X1(a,b,j,d) --> scr1(a,b,d,j)
c--------------------------------------------------------------------------

               next = 0
               do g = g1, g2
               do h = h1, h2
               do f = f1, f2
               do e = e1, e2
                  next = next + 1 
                  scr1(next) = x1(e,f,g,h)
               enddo
               enddo
               enddo
               enddo
                
c---------------------------------------------------------------------------
c   Now, do scr1(a,b,d,j) * x2(c,j) --> scr2(a,b,d,c)
c---------------------------------------------------------------------------

               m = (e2-e1+1)*(f2-f1+1)*(h2-h1+1)
               n = (i2-i1+1)
               k = j2-j1+1
               lda = (e2-e1+1)*(f2-f1+1)*(h2-h1+1)
               ldb = (ni2-ni1+1)
               ldc = (a2-a1+1)*(b2-b1+1)*(d2-d1+1)
               call dgemm('N', 'T', m,n,k,1.d0, scr1, lda,
     *                        x2(i1,j1), ldb, 0.d0,
     *                        scr2, ldc)

c---------------------------------------------------------------------------
c   Now sort scr2(a,b,d,c) --> y(a,b,c,d)
c---------------------------------------------------------------------------

               next = 0
               do c = c1, c2
               do d = d1, d2
               do b = b1, b2
               do a = a1, a2
                  next = next + 1
                  y(a,b,c,d) = scr2(next)
               enddo
               enddo
               enddo
               enddo

               return
         endif

         if (icontract .eq. 4) then 

c--------------------------------------------------------------------------
c   Contract over the 4th index of X1.
c
c   Y(a,b,c,d) = X1(a,b,c,j) * X2(d,j)
c--------------------------------------------------------------------------

               m = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)
               n = d2-d1+1
               k = j2-j1+1
               lda = (ne2-ne1+1)*(nf2-nf1+1)*(ng2-ng1+1)
               ldb = ni2-ni1+1
               ldc = (na2-na1+1)*(nb2-nb1+1)*(nc2-nc1+1)
               call dgemm('N', 'T', m,n,k,1.d0, x1(a1,b1,c1,j1),
     *                 lda, x2(i1,j1), ldb,
     *                 0.d0, y(a1,b1,c1,d1), ldc)
               return
         endif

      endif

      return
      end
