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
      subroutine tensor_contraction(x1handle, x1_block, x1_ndx,
     *                      x2handle, x2_block, x2_ndx,
     *                      yhandle, y_block, yndx,
     *                      array_table, narray_table, 
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table, 
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
c--------------------------------------------------------------------------
c   Forms the tensor contraction of blocks x1 and x2, resulting in block y.
c   x1, x2, and y are the array_table entries for their respective arrays.
c   x1_block, x2_block, and y_block are the actual block indices of the
c   blocks to be summed.
c
c   The "cind" argument is an array of indices over which the contraction
c   is performed.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table
      integer nscalar_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table) 
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      integer x1handle, x2handle, yhandle
      integer x1_block, x2_block, y_block 
      integer x1_ndx, x2_ndx, yndx
      integer flopcount

      integer n1, n2, ny
      integer nw, iseg, range1, range2
      integer*8 iy, ix1, ix2, daddr
      integer*8 get_block_index, get_index_from_base
      integer stack
      integer i, j, k
      integer ind(4,3), cindx1(2), cindx2(2)
      integer x1_type, x2_type, ytype
      double precision value

      double precision y(1), x1(1), x2(2)
#ifdef ALTIX
      pointer (ptry, y)
      pointer (ptrx1, x1)
      pointer (ptrx2, x2)
#else
      common y, x1, x2
#endif

      integer indy(4), segy(4), ya(4), yb(4)
      integer indx1(4), segx1(4), x1a(4), x1b(4)
      integer indx2(4), segx2(4), x2a(4), x2b(4)
      integer numy(4), numx1(4), numx2(4)
      integer nya(4), nyb(4), nx1a(4), nx1b(4), nx2a(4), nx2b(4)

      flopcount = 0

#ifdef ALTIX
      ptry = dshptr
      ptrx1 = dshptr
      ptrx2 = dshptr
#endif

c------------------------------------------------------------------------
c   Chcek for any scalar operands.
c------------------------------------------------------------------------

      x1_type = array_table(c_array_type,x1handle)
      x2_type = array_table(c_array_type,x2handle)
      ytype = array_table(c_array_type,yhandle)
      if (x1_type .eq. scalar_value .or.
     *    x2_type .eq. scalar_value .or.
     *    ytype .eq. scalar_value) then
         print *,'Error: tensor_contraction was called with ',
     *          'one or more scalar operands.'
         print *,'Operands are ',x1handle, x2handle, yhandle 
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Block contraction operation.
c-------------------------------------------------------------------------

      n1 = array_table(c_nindex, x1handle)
      n2 = array_table(c_nindex, x2handle)
      ny = array_table(c_nindex, yhandle)

c---------------------------------------------------------------------------
c   Get the array segments and ranges of each block.
c---------------------------------------------------------------------------

      do i = 1, ny
         indy(i) = array_table(c_index_array1+i-1, yhandle)
         segy(i) = index_table(c_current_seg,indy(i))
         call get_index_segment(indy(i), segy(i), segment_table,
     *                     nsegment_table, index_table,
     *                     nindex_table, ya(i), yb(i))  
         numy(i) = yb(i) - ya(i) + 1
      enddo

      stack = array_table(c_array_stack,yhandle)
      iy = get_block_index(yhandle, y_block, stack, yndx, y, .true.)
      if (iy .eq. -1) then
         daddr = address_table(yhandle)
         if (daddr .eq. 0) then
            print *,'Error: Cannot form address for array ',yhandle,
     *                ' block ',y_block,'.'
            call dump_block_ids()
            call abort_job()
         endif

         iy = get_index_from_base(daddr, y, 2)

         do i = 1, ny
            nya(i) = array_table(c_index_range1+i-1,yhandle)
            nyb(i) = array_table(c_index_range2+i-1,yhandle)
         enddo
      else
         do i = 1, ny
            nya(i) = ya(i)
            nyb(i) = yb(i)
         enddo
      endif

      do i = 1, n1
         indx1(i) = array_table(c_index_array1+i-1, x1handle)
         segx1(i) = index_table(c_current_seg,indx1(i))
         call get_index_segment(indx1(i), segx1(i), segment_table,
     *                     nsegment_table, index_table,
     *                     nindex_table, x1a(i), x1b(i))  
         numx1(i) = x1b(i) - x1a(i) + 1
      enddo

      stack = array_table(c_array_stack,x1handle)
      ix1 = get_block_index(x1handle, x1_block, stack, 
     *                      x1_ndx, x1, .true.)
      if (ix1 .eq. -1) then
         daddr = address_table(x1handle)
         if (daddr .eq. 0) then
            print *,'Error: Cannot form address for array ',x1handle,
     *                ' block ',x1_block,'.'
            call dump_block_ids()
            call abort_job()
         endif

         ix1 = get_index_from_base(daddr, x1, 2)
         do i = 1, n1
            nx1a(i) = array_table(c_index_range1+i-1,x1handle)
            nx1b(i) = array_table(c_index_range2+i-1,x1handle)
         enddo
      else
         do i = 1, n1
            nx1a(i) = x1a(i)
            nx1b(i) = x1b(i)
         enddo
      endif

      do i = 1, n2
         indx2(i) = array_table(c_index_array1+i-1, x2handle)
         segx2(i) = index_table(c_current_seg,indx2(i))
         call get_index_segment(indx2(i), segx2(i), segment_table,
     *                     nsegment_table, index_table,
     *                     nindex_table, x2a(i), x2b(i))  
         numx2(i) = x2b(i) - x2a(i) + 1
      enddo

      stack = array_table(c_array_stack,x2handle)
      ix2 = get_block_index(x2handle, x2_block, stack,
     *                      x2_ndx, x2, .true.)
      if (ix2 .eq. -1) then
         daddr = address_table(x2handle)
         if (daddr .eq. 0) then
            print *,'Error: Cannot form address for array ',x2handle,
     *                ' block ',x2_block,'.'
            call dump_block_ids()
            call abort_job()
         endif

         ix2 = get_index_from_base(daddr, x2, 2)
         do i = 1, n2
            nx2a(i) = array_table(c_index_range1+i-1,x2handle)
            nx2b(i) = array_table(c_index_range2+i-1,x2handle)

         enddo
      else
         do i = 1, n2
            nx2a(i) = x2a(i)
            nx2b(i) = x2b(i)
         enddo
      endif

      if (ny .ne. 4 .or. n1 .ne. 2 .or. n2 .ne. 2) then
         print *,'Error: Improper array dimensions for ',
     *            ' tensor_contraction: got ny, n1, n2 = ',
     *            ny,n1,n2,' should be 4, 2, 2.'
         print *,'Array handles: ',yhandle, x1handle, x2handle
         print *,' y: ',yhandle,y_block, yndx
         print *,' x1: ',x1handle, x1_block,x1_ndx
         print *,' x2: ',x2handle, x2_block,x2_ndx
         call abort_job()
      endif

      call twork4222(y(iy),nya(1),nyb(1),nya(2),nyb(2),
     *                     nya(3),nyb(3),nya(4),nyb(4),
     *                     ya(1),yb(1),ya(2),yb(2),ya(3),yb(3),
     *                     ya(4),yb(4),indy,
     *                     x1(ix1),nx1a(1),nx1b(1), nx1a(2),nx1b(2),
     *                     x1a(1),x1b(1),x1a(2),x1b(2), indx1,
     *                     x2(ix2),nx2a(1),nx2b(1), nx2a(2),nx2b(2),
     *                     x2a(1),x2b(1),x2a(2),x2b(2),indx2,
     *                     flopcount)

      return
      end
