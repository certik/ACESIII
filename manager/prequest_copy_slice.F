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
      subroutine prequest_copy_slice(node, source, target,
     *                server_table, nserver_table)
c-------------------------------------------------------------------------
c   Copies a slice of data from "source" to "target" based on description
c   of data contained in a message node of the server.
c-------------------------------------------------------------------------
      implicit none
      include 'server.h'

      integer i, node
      integer nserver_table
      integer server_table(*)
      double precision source(*), target(*)
 
      integer nind
      integer bsrc(mx_array_index), esrc(mx_array_index)
      integer btar(mx_array_index), etar(mx_array_index)

c-------------------------------------------------------------------------
c   Determine the number of indices.
c-------------------------------------------------------------------------

      nind = server_msg(c_msg_nind,node)
    
c--------------------------------------------------------------------------
c   Gather the source and target segment ranges.
c--------------------------------------------------------------------------

      do i = 1, nind
         bsrc(i) = server_msg(c_msg_bsegs2+i-1,node) 
         esrc(i) = server_msg(c_msg_esegs2+i-1,node) 
         btar(i) = server_msg(c_msg_bsegs+i-1,node) 
         etar(i) = server_msg(c_msg_esegs+i-1,node) 
      enddo
      
c---------------------------------------------------------------------------
c   Call the appropriate routine to handle the copy.
c---------------------------------------------------------------------------

      if (nind .eq. 1) then
         call preq_copy1(source, bsrc(1), esrc(1),
     *                   target, btar(1), etar(1))  
      else if (nind .eq. 2) then
         call preq_copy2(source, bsrc(1), esrc(1), bsrc(2), esrc(2),
     *                   target, btar(1), etar(1), btar(2), etar(2))
      else if (nind .eq. 3) then
         call preq_copy3(source, bsrc(1), esrc(1), bsrc(2), esrc(2),
     *                           bsrc(3), esrc(3),
     *                   target, btar(1), etar(1), btar(2), etar(2),
     *                           btar(3), etar(3))
      else if (nind .eq. 4) then
         call preq_copy4(source, bsrc(1), esrc(1), bsrc(2), esrc(2),
     *                           bsrc(3), esrc(3), bsrc(4), esrc(4),
     *                   target, btar(1), etar(1), btar(2), etar(2),
     *                           btar(3), etar(3),  btar(4), etar(4))
      else if (nind .eq. 5) then
         call preq_copy5(source, bsrc(1), esrc(1), bsrc(2), esrc(2),
     *                           bsrc(3), esrc(3), bsrc(4), esrc(4),
     *                           bsrc(5), esrc(5),
     *                   target, btar(1), etar(1), btar(2), etar(2),
     *                           btar(3), etar(3), btar(4), etar(4),
     *                           btar(5), etar(5))
      else if (nind .eq. 6) then
         call preq_copy6(source, bsrc(1), esrc(1), bsrc(2), esrc(2),
     *                           bsrc(3), esrc(3), bsrc(4), esrc(4),
     *                           bsrc(5), esrc(5), bsrc(6), esrc(6),
     *                   target, btar(1), etar(1), btar(2), etar(2),
     *                           btar(3), etar(3), btar(4), etar(4),
     *                           btar(5), etar(5), btar(6), etar(6))
      else
         print *,'Error: PREQUEST can only handle up to 6 indices.'
         call server_abort_job(server_table, nserver_table)
      endif
      return
      end

      subroutine preq_copy1(x, a1, a2, y, i1, i2)
      implicit none
      integer a1, a2, i1, i2
      double precision x(a1:a2)
      double precision y(i1:i2)

      integer a, i

      do i = i1, i2
         y(i) = x(i)
      enddo
     
      return
      end

      subroutine preq_copy2(x, a1, a2, b1, b2, 
     *                      y, i1, i2, j1, j2)
      implicit none
      integer a1, a2, b1, b2,
     *        i1, i2, j1, j2
      double precision x(a1:a2,b1:b2)
      double precision y(i1:i2,j1:j2)

      integer a, b, i, j

      do j = j1, j2
      do i = i1, i2
         y(i,j) = x(i,j)
      enddo
      enddo

      return
      end

      subroutine preq_copy3(x, a1, a2, b1, b2, c1, c2,
     *                      y, i1, i2, j1, j2, k1, k2)
      implicit none
      integer a1, a2, b1, b2, c1, c2,
     *        i1, i2, j1, j2, k1, k2
      double precision x(a1:a2,b1:b2,c1:c2)
      double precision y(i1:i2,j1:j2,k1:k2)

      integer a, b, c, i, j, k

      do k = k1, k2
      do j = j1, j2
      do i = i1, i2
         y(i,j,k) = x(i,j,k)
      enddo
      enddo
      enddo

      return
      end

      subroutine preq_copy4(x, a1, a2, b1, b2, c1, c2, d1, d2,
     *                      y, i1, i2, j1, j2, k1, k2, l1, l2)
      implicit none
      integer a1, a2, b1, b2, c1, c2, d1, d2,
     *        i1, i2, j1, j2, k1, k2, l1, l2
      double precision x(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision y(i1:i2,j1:j2,k1:k2,l1:l2)

      integer a, b, c, d, i, j, k, l

      do l = l1, l2
      do k = k1, k2
      do j = j1, j2
      do i = i1, i2
         y(i,j,k,l) = x(i,j,k,l)
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine preq_copy5(x, a1, a2, b1, b2, c1, c2, d1, d2, e1, e2,
     *                      y, i1, i2, j1, j2, k1, k2, l1, l2, m1, m2)
      implicit none
      integer a1, a2, b1, b2, c1, c2, d1, d2, e1, e2,
     *        i1, i2, j1, j2, k1, k2, l1, l2, m1, m2
      double precision x(a1:a2,b1:b2,c1:c2,d1:d2,e1:e2)
      double precision y(i1:i2,j1:j2,k1:k2,l1:l2,m1:m2)

      integer a, b, c, d, e, i, j, k, l, m

      do m = m1, m2
      do l = l1, l2
      do k = k1, k2
      do j = j1, j2
      do i = i1, i2
         y(i,j,k,l,m) = x(i,j,k,l,m)
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine preq_copy6(x, a1, a2, b1, b2, c1, c2, d1, d2, e1, e2,
     *                         f1, f2,
     *                      y, i1, i2, j1, j2, k1, k2, l1, l2, m1, m2,
     *                         n1, n2)
      implicit none
      integer a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2,
     *        i1, i2, j1, j2, k1, k2, l1, l2, m1, m2, n1, n2
      double precision x(a1:a2,b1:b2,c1:c2,d1:d2,e1:e2,f1:f2)
      double precision y(i1:i2,j1:j2,k1:k2,l1:l2,m1:m2,n1:n2)

      integer a, b, c, d, e, f, i, j, k, l, m, n

      do n = n1, n2
      do m = m1, m2
      do l = l1, l2
      do k = k1, k2
      do j = j1, j2
      do i = i1, i2
         y(i,j,k,l,m,n) = x(i,j,k,l,m,n)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

