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
      subroutine cwork1322(y,na1,na2,a1,a2,inda,
     *                     x1,nd1,nd2,ne1,ne2,nf1,nf2,
     *                     d1,d2,e1,e2,f1,f2, indx1,
     *                     x2,ni1,ni2,nj1,nj2,
     *                     i1,i2,j1,j2,indi, indj,
     *                     scr1, scr2,flopcount) 
c-------------------------------------------------------------------------
c   Performs a "1322" contraction: 
c      1 index output array  --> y 
c      3 index operand array --> x1 
c      2 index operand array --> x2 
c      2 index contraction.
c
c   I. e., contract the 2 indices of the 2nd operand array out of the 
c   first operand array into the 2-index output array.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'

      character*1 char(5)  
      character*2 char2a, char2b   
      character*5 char5 

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,i1,i2,j1,j2
      integer inda, indb, indc, indx1(3), indi, indj, h(3), perm(5) 
      integer flopcount

      double precision y(na1:na2)
      double precision x1(nd1:nd2,ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2)
      double precision scr1(*), scr2(*), x  

      integer ia, ib, ii, ij, ix(3)
      integer i, j, a, b, c, d, e, f 
#ifdef __crayx1
      integer*8 ncount, m, n, k, lda, ldb, ldc 
#else
      integer ncount, m, n, k, lda, ldb, ldc 
#endif

      logical check
      check = .false.

c determine first uncontracted index

      do i = 1, 3 
         if (inda .eq. indx1(i)) perm(i) = 1 
      enddo 

      if (inda .eq. indi) perm(4) = 1 
      if (inda .eq. indj) perm(5) = 1 

      write(6,*) ' Entering cwork1322 ' 

c determine contracted index 

      do i = 1, 3 
         if (indx1(i) .eq. indi) then 
            perm(i) = 4 
            perm(4) = 4 
            go to 33 
         endif 
      enddo 

33    continue 

      do i = 1, 3 
         if (indx1(i) .eq. indj) then 
            perm(i) = 5 
            perm(5) = 5 
            go to 44 
         endif 
      enddo 

44    continue 

c classify permutation 

      char(1) = '1'
      char(2) = '2'
      char(3) = '3'
      char(4) = '4'
      char(5) = '5'

      char5 = char(perm(1))//char(perm(2))//char(perm(3))//char(perm(4))
     *                     //char(perm(5)) 

c consider all possible combinations in computation 

c consider all possible combinations in computation

      if (char5 .eq. '14545') then
c                     aijij
c
         do a = a1, a2
            x = 0.0d0
            do i = i1, i2
            do j = j1, j2
               x = x + x1(a,i,j)*x2(i,j)
            enddo
            enddo
            y(a) = x
         enddo

         go to 66
c
      endif

      if (char5 .eq. '45145') then
c                     ijaij
c
         do a = a1, a2
            x = 0.0d0
            do i = i1, i2
            do j = j1, j2
               x = x + x1(i,j,a)*x2(i,j)
            enddo
            enddo
            y(a) = x
         enddo

         go to 66
c
      endif


c
      write(6,*) ' THE PERMUTATION ', char5, ' not implemented in 
     *             cwork1322' 
       call abort_job()
66    continue 

      do a = a1, a2 
         write(6,*) a, y(a) 
      enddo 

      write(6,*) ' Leaving cwork1322 ' 

      return
      end

