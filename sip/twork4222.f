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
      subroutine twork4222(y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                     a1,a2,b1,b2,c1,c2,d1,d2,indy,
     *                     x1,ne1,ne2,nf1,nf2,
     *                     e1,e2,f1,f2, indx1,
     *                     x2,ni1,ni2,nj1,nj2,
     *                     i1,i2,j1,j2,indx2,
     *                     flopcount) 
c-------------------------------------------------------------------------
c   Performs a "4222" tensor contraction: 
c      4 index output array 
c      2 index operand array
c      2 index operand array
c      2 index contraction.
c
c--------------------------------------------------------------------------
      implicit none

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,nf1,nf2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,i1,i2,j1,j2
      integer indy(4), indx1(2), indx2(2)
      integer flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision x1(ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2)

      integer ia, ib, ic, id, ia2, ib2, ic2, id2
      integer ix(0:2), ix2(0:2)
      integer i, j, a, b, c, d

c---------------------------------------------------------------------------
c   Find which indices of the "x1" operand match the various y
c   indices.
c---------------------------------------------------------------------------

      do i = 0,2 
         ix(i) = 0
         ix2(i) = 0
      enddo

      ia = 0
      ib = 0
      ic = 0
      id = 0
      do j = 1, 2
      do i = 1, 4
         if (indx1(j) .eq. indy(i)) then
            if (i .eq. 1) then
               ia = indy(1)
            else if (i .eq. 2) then
               ib = indy(2)
            else if (i .eq. 3) then
               ic = indy(3)
            else 
               id = indy(4)
            endif
         endif 
      enddo
      enddo

      ia2 = 0
      ib2 = 0
      ic2 = 0
      id2 = 0
      do j = 1, 2
      do i = 1, 4
         if (indx2(j) .eq. indy(i)) then
            if (i .eq. 1) then
               ia2 = indy(1)
            else if (i .eq. 2) then
               ib2 = indy(2)
            else if (i .eq. 3) then
               ic2 = indy(3)
            else 
               id2 = indy(4)
            endif
         endif 
      enddo
      enddo

      flopcount = (a2-a1+1)*(b2-b1+1)*(c2-c1+1)*(d2-d1+1)
     
c----------------------------------------------------------------------------
c   ia is the 1st index of x1.
c----------------------------------------------------------------------------

      if (indx1(1) .eq. ia .and. indx1(2) .eq. ib .and.
     *    indx2(1) .eq. ic2 .and. indx2(2) .eq. id2) then
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d) = x1(a,b)*x2(c,d)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ia .and. indx1(2) .eq. ib .and.
     *    indx2(1) .eq. id2 .and. indx2(2) .eq. ic2) then
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d) = x1(a,b)*x2(d,c)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. ia .and. indx1(2) .eq. ic .and.
     *    indx2(1) .eq. ib2 .and. indx2(2) .eq. id2) then
         do d = d1,d2
         do b = b1,b2
         do c = c1,c2
         do a = a1,a2
            y(a,b,c,d) = x1(a,c)*x2(b,d)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ia .and. indx1(2) .eq. ic .and.
     *    indx2(1) .eq. id2 .and. indx2(2) .eq. ib2) then
         do d = d1,d2
         do b = b1,b2
         do c = c1,c2
         do a = a1,a2
            y(a,b,c,d) = x1(a,c)*x2(d,b)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. ia .and. indx1(2) .eq. id .and.
     *    indx2(1) .eq. ib2 .and. indx2(2) .eq. ic2) then
         do c = c1,c2
         do b = b1,b2
         do d = d1,d2
         do a = a1,a2
            y(a,b,c,d) = x1(a,d)*x2(b,c)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ia .and. indx1(2) .eq. id .and.
     *    indx2(1) .eq. ic2 .and. indx2(2) .eq. ib2) then
         do c = c1,c2
         do b = b1,b2
         do d = d1,d2
         do a = a1,a2
            y(a,b,c,d) = x1(a,d)*x2(c,b)
         enddo
         enddo
         enddo
         enddo

         return
      endif

c-------------------------------------------------------------------------
c   ib is 1st index in x1.
c-------------------------------------------------------------------------

      if (indx1(1) .eq. ib .and. indx1(2) .eq. ia .and.
     *    indx2(1) .eq. ic2 .and. indx2(2) .eq. id2) then
         do d = d1,d2
         do c = c1,c2
         do a = a1,a2
         do b = b1,b2
            y(a,b,c,d) = x1(b,a)*x2(c,d)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ib .and. indx1(2) .eq. ia .and.
     *    indx2(1) .eq. id2 .and. indx2(2) .eq. ic2) then
         do d = d1,d2
         do c = c1,c2
         do a = a1,a2
         do b = b1,b2
            y(a,b,c,d) = x1(b,a)*x2(d,c)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. ib .and. indx1(2) .eq. ic .and.
     *    indx2(1) .eq. ia2 .and. indx2(2) .eq. id2) then
         do d = d1,d2
         do a = a1,a2
         do c = c1,c2
         do b = b1,b2
            y(a,b,c,d) = x1(b,c)*x2(a,d)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ib .and. indx1(2) .eq. ic .and.
     *    indx2(1) .eq. id2 .and. indx2(2) .eq. ia2) then
         do a = a1,a2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
            y(a,b,c,d) = x1(b,c)*x2(d,a)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. ib .and. indx1(2) .eq. id .and.
     *    indx2(1) .eq. ia2 .and. indx2(2) .eq. ic2) then
         do c = c1,c2
         do a = a1,a2
         do d = d1,d2
         do b = b1,b2
            y(a,b,c,d) = x1(b,d)*x2(a,c)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ib .and. indx1(2) .eq. id .and.
     *    indx2(1) .eq. ic2 .and. indx2(2) .eq. ia2) then
         do a = a1,a2
         do c = c1,c2
         do d = d1,d2
         do b = b1,b2
            y(a,b,c,d) = x1(b,d)*x2(c,a)
         enddo
         enddo
         enddo
         enddo

         return
      endif

c----------------------------------------------------------------------------
c   ic is the 1st index of x1.
c----------------------------------------------------------------------------

      if (indx1(1) .eq. ic .and. indx1(2) .eq. ib .and.
     *    indx2(1) .eq. ia2 .and. indx2(2) .eq. id2) then
         do d = d1,d2
         do a = a1,a2
         do b = b1,b2
         do c = c1,c2
            y(a,b,c,d) = x1(c,b)*x2(a,d)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ic .and. indx1(2) .eq. ib .and.
     *    indx2(1) .eq. id2 .and. indx2(2) .eq. ia2) then
         do a = a1,a2
         do d = d1,d2
         do b = b1,b2
         do c = c1,c2
            y(a,b,c,d) = x1(c,b)*x2(d,a)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. ic .and. indx1(2) .eq. ia .and.
     *    indx2(1) .eq. ib2 .and. indx2(2) .eq. id2) then
         do d = d1,d2
         do b = b1,b2
         do a = a1,a2
         do c = c1,c2
            y(a,b,c,d) = x1(c,a)*x2(b,d)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ic .and. indx1(2) .eq. ia .and.
     *    indx2(1) .eq. id2 .and. indx2(2) .eq. ib2) then
         do d = d1,d2
         do b = b1,b2
         do a = a1,a2
         do c = c1,c2
            y(a,b,c,d) = x1(c,a)*x2(d,b)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. ic .and. indx1(2) .eq. id .and.
     *    indx2(1) .eq. ib2 .and. indx2(2) .eq. ia2) then
         do a = a1,a2
         do b = b1,b2
         do d = d1,d2
         do c = c1,c2
            y(a,b,c,d) = x1(c,d)*x2(b,a)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. ic .and. indx1(2) .eq. id .and.
     *    indx2(1) .eq. ia2 .and. indx2(2) .eq. ib2) then
         do b = b1,b2
         do a = a1,a2
         do d = d1,d2
         do c = c1,c2
            y(a,b,c,d) = x1(c,d)*x2(a,b)
         enddo
         enddo
         enddo
         enddo

         return
      endif

c----------------------------------------------------------------------------
c   id is the 1st index of x1.
c----------------------------------------------------------------------------

      if (indx1(1) .eq. id .and. indx1(2) .eq. ib .and.
     *    indx2(1) .eq. ic2 .and. indx2(2) .eq. ia2) then
         do a = a1,a2
         do c = c1,c2
         do b = b1,b2
         do d = d1,d2
            y(a,b,c,d) = x1(d,b)*x2(c,a)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. id .and. indx1(2) .eq. ib .and.
     *    indx2(1) .eq. ia2 .and. indx2(2) .eq. ic2) then
         do c = c1,c2
         do a = a1,a2
         do b = b1,b2
         do d = d1,d2
            y(a,b,c,d) = x1(d,b)*x2(a,c)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. id .and. indx1(2) .eq. ic .and.
     *    indx2(1) .eq. ib2 .and. indx2(2) .eq. ia2) then
         do a = a1,a2
         do b = b1,b2
         do c = c1,c2
         do d = d1,d2
            y(a,b,c,d) = x1(d,c)*x2(b,a)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. id .and. indx1(2) .eq. ic .and.
     *    indx2(1) .eq. ia2 .and. indx2(2) .eq. ib2) then
         do b = b1,b2
         do a = a1,a2
         do c = c1,c2
         do d = d1,d2
            y(a,b,c,d) = x1(d,c)*x2(a,b)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      if (indx1(1) .eq. id .and. indx1(2) .eq. ia .and.
     *    indx2(1) .eq. ib2 .and. indx2(2) .eq. ic2) then
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
         do d = d1,d2
            y(a,b,c,d) = x1(d,a)*x2(b,c)
         enddo 
         enddo 
         enddo 
         enddo 

         return
      endif

      if (indx1(1) .eq. id .and. indx1(2) .eq. ia .and.
     *    indx2(1) .eq. ic2 .and. indx2(2) .eq. ib2) then
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
         do d = d1,d2
            y(a,b,c,d) = x1(d,a)*x2(c,b)
         enddo
         enddo
         enddo
         enddo

         return
      endif

      print *,'twork4222: Your tensor operation does not fit',
     *    ' any valid pattern'
      call abort_job() 

      return
      end

