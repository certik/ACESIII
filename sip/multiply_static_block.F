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
      subroutine multiply_static_block(x, xblk, xndx, value, 
     *                         y, yblk, yndx,
     *                         array_table, narray_table, 
     *                         index_table, nindex_table,
     *                         segment_table, nsegment_table,
     *                         address_table)
c--------------------------------------------------------------------------
c   Multiply a block by a scalar constant.  Either the block is static or
c   the output block is static.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer x, xblk, y, yblk,narray_table, nindex_table,nsegment_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer*8 address_table(narray_table)
      integer xndx, yndx
      double precision value

      integer i, ind, seg, nindex
      integer y1(mx_array_index), y2(mx_array_index)
      integer ny1(mx_array_index), ny2(mx_array_index)
      integer nx1(mx_array_index), nx2(mx_array_index)
      double precision out(1), in(1)
#ifdef ALTIX
      pointer (ptr_out, out)
      pointer (ptr_in, in)
#else
      common out, in
#endif

      integer*8 indy, indx, daddr
      integer*8 get_index_from_base
      integer*8 get_block_index
      integer stack

#ifdef ALTIX
      ptr_out = dshptr
      ptr_in  = dshptr
#endif

      nindex = array_table(c_nindex,x)
      do i = 1, nindex
         ind = array_table(c_index_array1+i-1,y)
         seg = index_table(c_current_seg,ind)
         call get_index_segment(ind, seg, segment_table,
     *                        nsegment_table, index_table,
     *                        nindex_table, y1(i), y2(i))
      enddo

      if (array_table(c_array_type, x) .eq. static_array) then
         daddr = address_table(x)
         indx = get_index_from_base(daddr, in, 2)
         do i = 1, nindex
            nx1(i) = array_table(c_index_range1+i-1,x)
            nx2(i) = array_table(c_index_range2+i-1,x)
         enddo
      else
         stack = array_table(c_array_stack,x) 
         indx = get_block_index(x, xblk, stack, xndx, in, .true.)
         do i = 1, nindex
            nx1(i) = y1(i)
            nx2(i) = y2(i)
         enddo
      endif

      if (array_table(c_array_type, y) .eq. static_array) then
         daddr = address_table(y)
         indy = get_index_from_base(daddr, out, 2)
         do i = 1, nindex
            ny1(i) = array_table(c_index_range1+i-1,y)
            ny2(i) = array_table(c_index_range2+i-1,y)
         enddo
      else
         stack = array_table(c_array_stack,y)
         indy = get_block_index(y, yblk, stack,yndx, out, .true.)
         do i = 1, nindex
            ny1(i) = y1(i)
            ny2(i) = y2(i)
         enddo
      endif

      if (nindex .eq. 4) then
         call multiply_block4(out(indy), y1, y2, ny1(1),ny2(1),
     *            ny1(2),ny2(2),ny1(3),ny2(3),ny1(4),ny2(4),
     *            in(indx), nx1(1),nx2(1),
     *            nx1(2),nx2(2),nx1(3),nx2(3),nx1(4),nx2(4),
     *            value)
      else
         call multiply_block2(out(indy), y1, y2, ny1(1),ny2(1),
     *            ny1(2),ny2(2), in(indx), nx1(1),nx2(1),
     *            nx1(2),nx2(2), value)
      endif

      return
      end

       
      subroutine multiply_block2(y, y1, y2, na1, na2, nb1, nb2, 
     *                           x, nc1, nc2, nd1, nd2,
     *                           value)
      implicit none
      integer y1(*), y2(*)
      integer na1, na2, nb1, nb2, nc1, nc2, nd1, nd2
      double precision y(na1:na2, nb1:nb2)
      double precision x(nc1:nc2, nd1:nd2)
      double precision value

      integer a, b

      do b = y1(2), y2(2)
      do a = y1(1), y2(1)
         y(a,b) = value*x(a,b)
      enddo
      enddo
      return
      end

      subroutine multiply_block4(y, y1, y2, na1, na2, nb1, nb2, 
     *                           nc1, nc2, nd1, nd2, x, ne1, ne2,
     *                           nf1, nf2, ng1, ng2, nh1, nh2,
     *                           value)
      implicit none
      integer y1(*), y2(*)
      integer na1, na2, nb1, nb2, nc1, nc2, nd1, nd2
      integer ne1, ne2, nf1, nf2, ng1, ng2, nh1, nh2
      double precision y(na1:na2, nb1:nb2, nc1:nc2, nd1:nd2)
      double precision x(ne1:ne2, nf1:nf2, ng1:ng2, nh1:nh2)
      double precision value

      integer a, b, c, d

      do d = y1(4), y2(4)
      do c = y1(3), y2(3)
      do b = y1(2), y2(2)
      do a = y1(1), y2(1)
         y(a,b,c,d) = value*x(a,b,c,d)
      enddo
      enddo
      enddo
      enddo
      return
      end

