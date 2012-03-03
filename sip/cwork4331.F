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
      subroutine cwork4331(y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2, 
     *                     a1,a2,b1,b2,c1,c2,d1,d2,inda, indb,
     *                     indc, indd,x1,ne1,ne2,nf1,nf2,ng1,ng2,
     *                     e1,e2,f1,f2,g1,g2,
     *                     indx1, x2,ni1,ni2,nj1,nj2,nk1,nk2,
     *                     i1,i2,j1,j2,k1,k2,indi,indj,indk,
     *                     cind, flopcount,scr1, scr2)
c-------------------------------------------------------------------------
c   Performs a "4331" contraction:
c      4 index output array
c      3 index operand array
c      3 index operand array
c      1 index contraction.
c
c   I. e., contract one of the the 2 indices of the 2nd operand array out
c   of the first operand array, replacing the index by the 2nd operand's
c   2nd index.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'
      include 'parallel_info.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2 
      integer ne1,ne2,nf1,nf2,ng1,ng2
      integer ni1, ni2, nj1, nj2, nk1, nk2 
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,g1,g2,
     *        i1,i2,j1,j2,k1,k2  
      integer inda, indb, indc, indd, indx1(3), indi, indj, indk  
      integer cind(2), flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision yS(na1:na2,nb1:nb2,nc1:nc2)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2)
      double precision x2S(ni1:ni2,nj1:nj2)

      double precision scr1(*)
      double precision scr2(*)
      integer next

      integer ia, ib, ic, id, icontract, ix(3), iy
      integer i, j, a, b, c, d  
      integer indy(4), map(4)
      integer ie, if, ig
      integer e, f, g 
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
      do i = 1, 3

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
            print *,'Error: Invalid index in cwork4331'
            print *,'X1 index is ',indx1(i),' y indices: ',inda, indb,
     *              indc,indd,    
     *              ' x2 indices: ',indi,indj,indk,' cind(1) = ',cind(1)
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

c     write(6,*) ' CIND = ', cind(1), 'indi = ', indi, 
c    *           'icontract = ', icontract 

c     write(6,*) ' indy(i):', (indy(i),i=1,4) 
      if (indi .eq. cind(1)) then

c--------------------------------------------------------------------------
c   Contract over the first index in X2.
c--------------------------------------------------------------------------

         flopcount = 2*(i2-i1+1)*(c2-c1+1)*(b2-b1+1)*(a2-a1+1)
         if (indy(1) .eq. 1 .and. indy(2) .eq. 2 .and.
     *       indy(3) .eq. 3) then !  .and. indy(4) .eq. 4) then  

            if (icontract .eq. 1) then

c--------------------------------------------------------------------------
c   Contract over the 1st index of X1.
c
c   Y(a,b,c,d) = X1(i,b,c)*X2(i,a,d)
c--------------------------------------------------------------------------

            do d = d1, d2 

               do i = i1, i2 
               do j = j1, j2 
                  x2S(i,j) = x2(i,j,d) 
               enddo 
               enddo 

               m = (a2-a1+1)
               n = (b2-b1+1)*(c2-c1+1)
               k = i2-i1+1
               lda = (ni2-ni1+1)
               ldb = ne2-ne1+1
               ldc = (na2-na1+1)
               call dgemm('T', 'N', m,n,k,1.d0, x2S(i1,j1), lda,
     *                 x1(i1,b1,c1), ldb,
     *                 0.d0, yS(a1,b1,c1), ldc)

               do a = a1, a2 
               do b = b1, b2 
               do c = c1, c2 
                  y(a,b,c,d) = yS(a,b,c)
               enddo 
               enddo 
               enddo 

            enddo ! d 

               return

            endif

            if (icontract .eq. 2) then
c              write(6,*) ' Entering dgemm --> icontract = 2' 

c--------------------------------------------------------------------------
c   Contract over the 2nd index of X1.
c
c   Y(a,b,c,d) = X1(a,i,c) * X2(i,b,d)
c--------------------------------------------------------------------------

            do d = d1, d2 

               do i = i1, i2 
               do j = j1, j2 
                  x2S(i,j) = x2(i,j,d) 
               enddo 
               enddo 

               m = (a2-a1+1)
               n = (b2-b1+1)
               k = i2-i1+1
               lda = (ne2-ne1+1)
               ldb = (ni2-ni1+1)
               ldc = (na2-na1+1)
               do c = c1, c2
                  call dgemm('N', 'N', m,n,k,1.d0, x1(a1,i1,c), lda,
     *                        x2S(i1,j1), ldb, 0.d0,
     *                        yS(a1,b1,c), ldc)
               enddo

               do a = a1, a2 
               do b = b1, b2 
               do c = c1, c2 
                  y(a,b,c,d) = yS(a,b,c)
               enddo 
               enddo 
               enddo 

            enddo ! d 

               return
            endif

            if (icontract .eq. 3) then

c--------------------------------------------------------------------------
c   Contract over the 3rd index of X1.
c
c   Y(a,b,c,d) + X1(a,b,i)*X2(i,c,d)
c--------------------------------------------------------------------------

            do d = d1, d2 

               do i = i1, i2 
               do j = j1, j2 
                  x2S(i,j) = x2(i,j,d) 
               enddo 
               enddo 

               m = (a2-a1+1)*(b2-b1+1)
               n = (c2-c1+1)
               k = i2-i1+1
               lda = (ne2-ne1+1)*(nf2-nf1+1)
               ldb = (ni2-ni1+1)
               ldc = (na2-na1+1)*(nb2-nb1+1)
                  call dgemm('N', 'N', m,n,k,1.d0, x1(a1,b1,i1), lda,
     *                        x2S(i1,j1), ldb, 0.d0,
     *                        yS(a1,b1,c1), ldc)

               do a = a1, a2 
               do b = b1, b2 
               do c = c1, c2 
                  y(a,b,c,d) = yS(a,b,c)
               enddo 
               enddo 
               enddo 

            enddo ! d 
            endif
            return

         endif 

      endif 

      print *,'Error: Invalid contraction at line ',current_line
      print *,'X1 index is ',(indx1(i),i=1,3),' y indices: ',inda, indb,
     *              indc,indd,    
     *              ' x2 indices: ',indi,indj,indk,' cind(1) = ',cind(1)
      call abort_job()  
111   continue 
      return
      end
       
