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
      subroutine assign_block(xhandle, x_block, xndx,
     *                      yhandle, y_block, yndx,
     *                      array_table, narray_table, 
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table, 
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
c--------------------------------------------------------------------------
c   Transfers the data elements of block x_block of array xhandle into 
c   block y_block of array yhandle.  The indices of both arrays must
c   match, but need not be in the same order.  If the indices are out of 
c   order, the data is transposed instead of simply moved.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nscalar_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      double precision scalar_table(nscalar_table)
      integer nsegment_table 
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer*8 address_table(narray_table)
      integer xhandle, yhandle
      integer x_block, y_block 
      integer xndx, yndx
      integer flopcount

      integer nx, ny
      integer nw, nwx, iseg 
      integer yrange1(mx_array_index), yrange2(mx_array_index)
      integer xrange1(mx_array_index), xrange2(mx_array_index)
      integer ydim1(mx_array_index), ydim2(mx_array_index)
      integer xdim1(mx_array_index), xdim2(mx_array_index)
      integer*8 indy, indx, daddr
      integer*8 get_block_data_index
      integer stack
      integer*8 get_index_from_base
      integer get_blkmgr_blocksize
      integer i, j, k
      integer xind(mx_array_index), yind(mx_array_index)
      integer blkmgr_blocksize
      integer xtype, ytype
      integer index, seg, nindex

      logical match, transpose, slice_insert
      logical ysub_lookup(mx_array_index), xsub_lookup(mx_array_index)

      double precision y(1), x(1)
#ifdef ALTIX
      pointer (xptr, x)
      pointer (yptr, y)
#else
      common y, x
#endif
      double precision xval

#ifdef ALTIX
      xptr = dshptr
      yptr = dshptr
#endif

      xtype = array_table(c_array_type, xhandle)
      ytype = array_table(c_array_type, yhandle)

      if (ytype .ne. scalar_value) then

c---------------------------------------------------------------------------
c   Check array indices for slice/insertion operation using subindices.
c   First, check the type of the array table index against that of the 
c   underlying data block.
c---------------------------------------------------------------------------

         nindex = array_table(c_nindex, yhandle)
         if (ytype .eq. static_array) then
            do i = 1, nindex
               yind(i) = array_table(c_index_original+i-1, yhandle)
            enddo
         else
            call get_block_indices(yndx, yind)
         endif

         slice_insert = .false.
         do i = 1, nindex
            index = array_table(c_index_array1+i-1,yhandle)
            if (index_table(c_index_type, index) .eq. subindex .and.
     *          index_table(c_index_type, yind(i)) .ne. subindex) then
               ysub_lookup(i) = .true.
               slice_insert = .true.
            else
               ysub_lookup(i) = .false. 
            endif 
         enddo

         if (xtype .ne. scalar_value) then
            if (xtype .eq. static_array) then
               do i = 1, nindex
                  xind(i) = array_table(c_index_original+i-1,xhandle)
               enddo 
            else
               call get_block_indices(xndx, xind)
            endif

            do i = 1, nindex
               index = array_table(c_index_array1+i-1,xhandle)
               if (index_table(c_index_type, index) .eq. subindex .and.
     *          index_table(c_index_type, xind(i)) .ne. subindex) then
                  xsub_lookup(i) = .true.
                  slice_insert = .true.  
               else
                  xsub_lookup(i) = .false.
               endif
            enddo
         endif   ! xtype .ne. scalar_value

c----------------------------------------------------------------------------
c   Slices and insertions are handled separately by subroutine 
c   slice_insert_block.
c----------------------------------------------------------------------------

         if (slice_insert) then
            call slice_insert_block(xhandle, x_block, xndx,
     *                      yhandle, y_block, yndx,
     *                      array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, yind, xind, 
     *                      ysub_lookup, xsub_lookup, ytype, xtype)
            return
         endif
      endif   ! ytype .ne. scalar_value

      if (xtype .eq. scalar_value) then
         i = array_table(c_scalar_index, xhandle)
         xval = scalar_table(i)

         if (ytype .eq. scalar_value) then
         
c--------------------------------------------------------------------------
c   Scalar-to-scalar assignment.
c--------------------------------------------------------------------------
 
            i = array_table(c_scalar_index,yhandle)
            scalar_table(i) = xval
            return
         else

c---------------------------------------------------------------------------
c   Scalar-to-block assignment.  Fill block with copies of xval.
c---------------------------------------------------------------------------

           if (ytype .eq. static_array) then
              daddr = address_table(yhandle)
              indy  = get_index_from_base(daddr, y, 2)
              nindex = array_table(c_nindex, yhandle)
              do i = 1, nindex
                 index = array_table(c_index_array1+i-1,yhandle)
                 seg   = index_table(c_current_seg, index)
                 call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, yrange1(i),
     *                             yrange2(i))
                 ydim1(i) = array_table(c_index_range1+i-1,yhandle)
                 ydim2(i) = array_table(c_index_range2+i-1,yhandle)
              enddo

              if (nindex .eq. 4) then
                 call fill_static_block4(y(indy),yrange1(1),yrange2(1), 
     *                     yrange1(2), yrange2(2),
     *                     yrange1(3), yrange2(3),
     *                     yrange1(4), yrange2(4),
     *                     ydim1(1),ydim2(1),ydim1(2),ydim2(2),
     *                     ydim1(3),ydim2(3),ydim1(4),ydim2(4),xval)
              else if (nindex .eq. 2) then
                 call fill_static_block2(y(indy),yrange1(1),yrange2(1),
     *                     yrange1(2), yrange2(2), 
     *                     ydim1(1),ydim2(1),ydim1(2),ydim2(2),
     *                     xval)
              else
                print *,'Error: assign scalar must have a target ',
     *                  'array of 2 or 4 indices.'
                call abort_job()    
                 
              endif
           else
              stack = array_table(c_array_stack,yhandle)
              call fill_block(yhandle, y_block, stack, yndx, xval)
           endif
           return
         endif 
      endif
 
      nx = array_table(c_nindex, xhandle)
      ny = array_table(c_nindex, yhandle)

      flopcount = 0

c-----------------------------------------------------------------------
c   Match the operand indices with that of the result.
c-----------------------------------------------------------------------

      do i = 1, ny
         yind(i) = array_table(c_index_array1+i-1,yhandle)
      enddo

      do i = 1, nx
         xind(i) = array_table(c_index_array1+i-1,xhandle)
      enddo

      if (ny .lt. nx) then ! VFL 
         
c--------------------------------------------------------------------------
c   Final indices of y must all be of type "index".
c--------------------------------------------------------------------------

         do i = ny+1, nx
            if (index_table(c_index_type, xind(i)) .ne. 
     *                                      simple_index) then
               print *,'Error: Final indices of destination array ',
     *                 'must all be of type index'
               write(6,*) ' Indices of Y:', (yind(k), k=1, ny) 
               write(6,*) ' Indices of X:', (xind(k), k=1, nx) 
               call abort_job() 
            endif
         enddo

         do i = 1, ny
            if (yind(i) .ne. xind(i)) then
               print *,'Error: First ',nx,' indices of arrays must ',
     *             ' match in assign statement'
               call abort_job()
            endif
         enddo
      else if (ny .gt. nx) then ! VFL 
         do i = nx+1, ny
            if (index_table(c_index_type, yind(i)) .ne.
     *                                      simple_index) then
               print *,'Error: Final indices of source array ',
     *                 'must all be of type index'
               write(6,*) ' Indices of Y:', (yind(k), k=1, ny) 
               write(6,*) ' Indices of X:', (xind(k), k=1, nx) 
               call abort_job()
            endif
         enddo

         do i = 1, nx
            if (yind(i) .ne. xind(i)) then
               print *,'Error: First ',ny,' indices of arrays must ',
     *             ' match in assign statement'
               call abort_job()
            endif
         enddo
      else   
         do i = 1, ny
            match = .false.
            do j = 1, nx
               if (yind(i) .eq. xind(j)) match = .true.
            enddo

            if (.not. match) then
               print *,'Error in assign_block: Non-matching indices ',
     *                  ' for arrays ',yhandle,' and ',xhandle 
               print *,'Array ',yhandle,' indices: ',
     *                   (yind(j),j=1,ny)
               print *,'Array ',xhandle,' indices: ',
     *                   (xind(j),j=1,nx)
               call abort_job()
            endif
         enddo
      endif

c--------------------------------------------------------------------------
c   Indices match ok, check for a straight copy or transpose.
c--------------------------------------------------------------------------

      transpose = .false.
      do i = 1, min(nx,ny)
         if (yind(i) .ne. xind(i) .and.
     *       .not. transpose) transpose = .true.
      enddo

      if (transpose .and. (nx .ne. ny)) then ! VFL 
         print *,'Transpose not allowed in assign_block if nx .ne. ny '
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Pick up the segment ranges, compute block size.
c--------------------------------------------------------------------------

      nw = 1
      do i = 1, ny
         iseg = index_table(c_current_seg,yind(i))
         call get_index_segment(yind(i), iseg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, yrange1(i), yrange2(i))
         nw = nw * (yrange2(i)-yrange1(i)+1)
      enddo

      nwx = 1 
      do i = 1, nx
         iseg = index_table(c_current_seg,xind(i))
         call get_index_segment(xind(i), iseg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, xrange1(i), xrange2(i))
         nwx = nwx * (xrange2(i)-xrange1(i)+1)
      enddo
      if (nw .ne. nwx) then 
         print *,'Number of elements being transfered in assign block
     *            incorrect', nw, nwx 
         call abort_job()
      endif

      if (ytype .eq. static_array) then
         do i = 1, ny
            ydim1(i) = array_table(c_index_range1+i-1,yhandle) 
            ydim2(i) = array_table(c_index_range2+i-1,yhandle)
         enddo

         daddr = address_table(yhandle)
         indy = get_index_from_base(daddr, y, 2) 
      else
         do i = 1, ny
            ydim1(i) = yrange1(i)
            ydim2(i) = yrange2(i)
         enddo
         stack = array_table(c_array_stack,yhandle)
         indy = get_block_data_index(yhandle, y_block, stack, yndx, y) 
      endif

      if (xtype .eq. static_array) then
         do i = 1, nx
            xdim1(i) = array_table(c_index_range1+i-1,xhandle)
            xdim2(i) = array_table(c_index_range2+i-1,xhandle)
         enddo

         daddr = address_table(xhandle)
         indx = get_index_from_base(daddr, x, 2) 
      else
         do i = 1, nx
            xdim1(i) = xrange1(i)
            xdim2(i) = xrange2(i)
         enddo
         stack = array_table(c_array_stack,xhandle)
         indx = get_block_data_index(xhandle, x_block, stack,
     *                               xndx, x) 
      endif

      if (transpose) then
         if (ny .eq. 4) then
            call transpose_4index(y(indy),
     *                   yrange1(1),yrange2(1), yrange1(2),yrange2(2),
     *                   yrange1(3),yrange2(3), yrange1(4),yrange2(4),
     *                   ydim1(1), ydim2(1),ydim1(2),ydim2(2),
     *                   ydim1(3), ydim2(3),ydim1(4),ydim2(4),
     *                x(indx),   
     *                   xrange1(1),xrange2(1), xrange1(2),xrange2(2),
     *                   xrange1(3),xrange2(3), xrange1(4),xrange2(4),
     *                   xdim1(1), xdim2(1),xdim1(2),xdim2(2),
     *                   xdim1(3), xdim2(3),xdim1(4),xdim2(4),
     *                 yind, xind)
      else if (ny .eq. 3) then
            call transpose_3index(y(indy),
     *                   yrange1(1),yrange2(1), yrange1(2),yrange2(2),
     *                   yrange1(3),yrange2(3), 
     *                   ydim1(1), ydim2(1),ydim1(2),ydim2(2),
     *                   ydim1(3), ydim2(3),
     *                x(indx),
     *                   xrange1(1),xrange2(1), xrange1(2),xrange2(2),
     *                   xrange1(3),xrange2(3), 
     *                   xdim1(1), xdim2(1),xdim1(2),xdim2(2),
     *                   xdim1(3), xdim2(3),
     *                 yind, xind)
      else if (ny .eq. 2) then

c------------------------------------------------------------------------
c   2 index array.
c------------------------------------------------------------------------

            call transpose_2index(y(indy),
     *                yrange1(1),yrange2(1), yrange1(2),yrange2(2),
     *                ydim1(1), ydim2(1),ydim1(2),ydim2(2),
     *                x(indx),
     *                xrange1(1),xrange2(1), xrange1(2),xrange2(2),
     *                xdim1(1), xdim2(1),xdim1(2),xdim2(2),
     *                yind, xind)
      else if (ny .eq. 6) then
            call transpose_6index(y(indy),
     *                   yrange1(1),yrange2(1), yrange1(2),yrange2(2),
     *                   yrange1(3),yrange2(3), yrange1(4),yrange2(4),
     *                   yrange1(5),yrange2(5), yrange1(6),yrange2(6),
     *                   ydim1(1), ydim2(1),ydim1(2),ydim2(2),
     *                   ydim1(3), ydim2(3),ydim1(4),ydim2(4),
     *                   ydim1(5), ydim2(5),ydim1(6),ydim2(6),
     *                x(indx),
     *                   xrange1(1),xrange2(1), xrange1(2),xrange2(2),
     *                   xrange1(3),xrange2(3), xrange1(4),xrange2(4),
     *                   xrange1(5),xrange2(5), xrange1(6),xrange2(6),
     *                   xdim1(1), xdim2(1),xdim1(2),xdim2(2),
     *                   xdim1(3), xdim2(3),xdim1(4),xdim2(4),
     *                   xdim1(5), xdim2(5),xdim1(6),xdim2(6),
     *                 yind, xind) 
         else
            print *,'Error: Incorrect number of indices for assign'
            call abort_job()
         endif
      else    ! .not. transpose
         if (ny .eq. 2) then
            call transpose_2index(y(indy),
     *                yrange1(1),yrange2(1), yrange1(2),yrange2(2),
     *                ydim1(1), ydim2(1),ydim1(2),ydim2(2),
     *                x(indx),
     *                xrange1(1),xrange2(1), xrange1(2),xrange2(2),
     *                xdim1(1), xdim2(1),xdim1(2),xdim2(2),
     *                yind, xind)
         else

c---------------------------------------------------------------------------
c   Straight copy.  Move the block.
c---------------------------------------------------------------------------

            do i = 1, nw
               y(indy+i-1) = x(indx+i-1)
            enddo
         endif
      endif

      return
      end

      subroutine transpose_4index(
     *                y, ya1,ya2,yb1,yb2,yc1,yc2,yd1,yd2,
     *                nya1,nya2,nyb1,nyb2,nyc1,nyc2,nyd1,nyd2,  
     *                x, xa1,xa2,xb1,xb2,xc1,xc2,xd1,xd2,
     *                nxa1,nxa2,nxb1,nxb2,nxc1,nxc2,nxd1,nxd2, 
     *                yind, xind) 
c--------------------------------------------------------------------------
c   Transposes data into from one array into another, using the block 
c   ranges in yrange and xrange.
c--------------------------------------------------------------------------
      implicit none

      integer ya1,ya2,yb1,yb2,yc1,yc2,yd1,yd2
      integer nya1,nya2,nyb1,nyb2,nyc1,nyc2,nyd1,nyd2
      integer xa1,xa2,xb1,xb2,xc1,xc2,xd1,xd2
      integer nxa1,nxa2,nxb1,nxb2,nxc1,nxc2,nxd1,nxd2
      double precision y(nya1:nya2,nyb1:nyb2,nyc1:nyc2,nyd1:nyd2)
      double precision x(nxa1:nxa2,nxb1:nxb2,nxc1:nxc2,nxd1:nxd2)
      integer yind(4), xind(4)

      integer i, j, n
      integer a,b,c,d
      integer xval(4), map(4)

      do i = 1, 4
         do j = 1, 4
            if (yind(i) .eq. xind(j)) map(i) = j
         enddo
      enddo

c--------------------------------------------------------------------------
c   Optimized transposes.
c--------------------------------------------------------------------------

      if (map(1) .eq. 1 .and. map(2) .eq. 2 .and. 
     *    map(3) .eq. 4 .and. map(4) .eq. 3) then
         n = (ya2-ya1+1)*(yb2-yb1+1)
         do c = yc1, yc2
         do d = yd1, yd2
            do b = yb1, yb2
            do a = ya1, ya2
               y(a,b,c,d) = x(a,b,d,c)
            enddo
            enddo
         enddo
         enddo
         return
      endif

      if (map(1) .eq. 2 .and. map(2) .eq. 1 .and.
     *    map(3) .eq. 3 .and. map(4) .eq. 4) then
         do d = yd1, yd2
         do c = yc1, yc2
         do a = ya1, ya2
            do b = yb1, yb2
               y(a,b,c,d) = x(b,a,c,d)
            enddo
         enddo
         enddo
         enddo

         return
      endif

c--------------------------------------------------------------------------
c   Generic 4-dimensional transpose.
c--------------------------------------------------------------------------

10101 continue
      do d = yd1, yd2
         xval(map(4)) = d 
      do c = yc1, yc2
         xval(map(3)) = c
      do b = yb1, yb2
         xval(map(2)) = b
      do a = ya1, ya2
         xval(map(1)) = a
         y(a,b,c,d) = x(xval(1), xval(2), xval(3), xval(4))
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine transpose_3index(
     *                y, ya1,ya2,yb1,yb2,yc1,yc2,
     *                nya1,nya2,nyb1,nyb2,nyc1,nyc2,  
     *                x, xa1,xa2,xb1,xb2,xc1,xc2,
     *                nxa1,nxa2,nxb1,nxb2,nxc1,nxc2,
     *                yind, xind) 
c--------------------------------------------------------------------------
c   Transposes data into from one array into another, using the block 
c   ranges in yrange and xrange.
c--------------------------------------------------------------------------
      implicit none

      integer ya1,ya2,yb1,yb2,yc1,yc2
      integer nya1,nya2,nyb1,nyb2,nyc1,nyc2
      integer xa1,xa2,xb1,xb2,xc1,xc2
      integer nxa1,nxa2,nxb1,nxb2,nxc1,nxc2
      double precision y(nya1:nya2,nyb1:nyb2,nyc1:nyc2)
      double precision x(nxa1:nxa2,nxb1:nxb2,nxc1:nxc2)
      integer yind(3), xind(3)

      integer i, j, n
      integer a,b,c,d
      integer xval(3), map(3)

      do i = 1, 3
         do j = 1, 3
            if (yind(i) .eq. xind(j)) map(i) = j
         enddo
      enddo

c--------------------------------------------------------------------------
c   Optimized transposes.
c--------------------------------------------------------------------------

10101 continue
      do c = yc1, yc2
         xval(map(3)) = c
      do b = yb1, yb2
         xval(map(2)) = b
      do a = ya1, ya2
         xval(map(1)) = a
         y(a,b,c) = x(xval(1), xval(2), xval(3))
      enddo
      enddo
      enddo

      return
      end

      subroutine transpose_6index(
     *                y, ya1,ya2,yb1,yb2,yc1,yc2,yd1,yd2,
     *                   ye1,ye2,yf1,yf2,
     *                nya1,nya2,nyb1,nyb2,nyc1,nyc2,nyd1,nyd2,  
     *                nye1,nye2,nyf1,nyf2,
     *                x, xa1,xa2,xb1,xb2,xc1,xc2,xd1,xd2,
     *                   xe1,xe2,xf1,xf2,
     *                nxa1,nxa2,nxb1,nxb2,nxc1,nxc2,nxd1,nxd2, 
     *                nxe1,nxe2,nxf1,nxf2,
     *                yind, xind) 
c--------------------------------------------------------------------------
c   Transposes data into from one array into another, using the block 
c   ranges in yrange and xrange.
c--------------------------------------------------------------------------
      implicit none

      integer ya1,ya2,yb1,yb2,yc1,yc2,yd1,yd2,ye1,ye2,yf1,yf2
      integer nya1,nya2,nyb1,nyb2,nyc1,nyc2,nyd1,nyd2
      integer nye1,nye2,nyf1,nyf2
      integer xa1,xa2,xb1,xb2,xc1,xc2,xd1,xd2,xe1,xe2,xf1,xf2
      integer nxa1,nxa2,nxb1,nxb2,nxc1,nxc2,nxd1,nxd2
      integer nxe1,nxe2,nxf1,nxf2
      double precision 
     *   y(nya1:nya2,nyb1:nyb2,nyc1:nyc2,nyd1:nyd2,nye1:nye2,nyf1:nyf2)
      double precision 
     *   x(nxa1:nxa2,nxb1:nxb2,nxc1:nxc2,nxd1:nxd2,nxe1:nxe2,nxf1:nxf2)
      integer yind(6), xind(6)

      integer i, j, n
      integer a,b,c,d,e,f
      integer xval(6), map(6)

      do i = 1, 6
         do j = 1, 6
            if (yind(i) .eq. xind(j)) map(i) = j
         enddo
      enddo

c--------------------------------------------------------------------------
c   Generic 4-dimensional transpose.
c--------------------------------------------------------------------------

10101 continue
      do f = yf1, yf2
         xval(map(6)) = f
      do e = ye1, ye2
         xval(map(5)) = e
      do d = yd1, yd2
         xval(map(4)) = d 
      do c = yc1, yc2
         xval(map(3)) = c
      do b = yb1, yb2
         xval(map(2)) = b
      do a = ya1, ya2
         xval(map(1)) = a
         y(a,b,c,d,e,f) = x(xval(1), xval(2), xval(3), xval(4), 
     *                      xval(5), xval(6))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine transpose_2index(
     *                y, ya1,ya2,yb1,yb2,
     *                nya1,nya2,nyb1,nyb2,
     *                x, xa1,xa2,xb1,xb2,
     *                nxa1,nxa2,nxb1,nxb2,
     *                yind, xind) 
c--------------------------------------------------------------------------
c   Transposes data into from one array into another, using the block 
c   ranges in yrange and xrange.
c--------------------------------------------------------------------------
      implicit none

      integer ya1,ya2,yb1,yb2
      integer nya1,nya2,nyb1,nyb2
      integer xa1,xa2,xb1,xb2
      integer nxa1,nxa2,nxb1,nxb2
      double precision y(nya1:nya2,nyb1:nyb2)
      double precision x(nxa1:nxa2,nxb1:nxb2)
      integer yind(2), xind(2)

      integer i, j
      integer a,b,c,d
      integer xval(2), map(2)

      do i = 1, 2
         do j = 1, 2
            if (yind(i) .eq. xind(j)) map(i) = j
         enddo
      enddo

      if (map(1) .eq. 1 .and. map(2) .eq. 2) then
         do b = yb1, yb2
         do a = ya1, ya2
            y(a,b) = x(a,b)
         enddo
         enddo
      else
         do a = ya1, ya2
         do b = yb1, yb2
            y(a,b) = x(b,a)
         enddo
         enddo

      endif

      return
      end

      subroutine fill_static_block4(y, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                             xval)
      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2
      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision xval

      integer a,b,c,d

      do d = d1,d2
      do c = c1,c2
      do b = b1,b2
      do a = a1,a2
         y(a,b,c,d) = xval
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine fill_static_block2(y, a1,a2,b1,b2,
     *                             na1,na2,nb1,nb2,
     *                             xval)
      integer a1,a2,b1,b2
      integer na1,na2,nb1,nb2
      double precision y(na1:na2,nb1:nb2)
      double precision xval

      integer a,b

      do b = b1,b2
      do a = a1,a2
         y(a,b) = xval
      enddo
      enddo

      return
      end

