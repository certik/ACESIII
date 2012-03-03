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
      subroutine slice_insert_block(xhandle, x_block, xndx,
     *                      yhandle, y_block, yndx,
     *                      array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, yind, xind,
     *                      ysub_lookup, xsub_lookup,
     *                      ytype, xtype)
c---------------------------------------------------------------------------
c   Handles subindex slice/insert operations from block x_block of array 
c   xhandle to block y_block of array yhandle.  The arrays yind and xind
c   have already been populated with the indices of their respective data
c   blocks, and arrays ysub_lookup and xsub_lookup are indicators telling
c   whether a subindex lookup for that index is required.
c----------------------------------------------------------------------------
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
      integer xind(mx_array_index), yind(mx_array_index)
      integer xtype, ytype
      
      logical ysub_lookup(mx_array_index), xsub_lookup(mx_array_index) 

      integer i, index, nindex, stack, seg, superindex
      integer*8 indy, indx, daddr
      integer*8 get_index_from_base
      integer*8 get_block_data_index
      integer ydim1(mx_array_index), ydim2(mx_array_index)
      integer xdim1(mx_array_index), xdim2(mx_array_index)
      integer yrange1(mx_array_index), yrange2(mx_array_index)
      integer xrange1(mx_array_index), xrange2(mx_array_index)
      integer blkseg(mx_array_index)
      integer ynindex, xnindex
      logical xscalar

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

c------------------------------------------------------------------------
c   Check for index mismatches.  One array may have more dimensions than 
c   the other, as long as all trailing indices are of type simple_index.
c------------------------------------------------------------------------

      xnindex = array_table(c_nindex, xhandle)
      ynindex = array_table(c_nindex, yhandle)
      if (xnindex .gt. ynindex) then
         nindex = ynindex
         if (xtype .ne. scalar_value) then
            do i = nindex+1, xnindex
               index = array_table(c_index_array1+i-1,xhandle)
               if (index_table(c_index_type,index) .ne. 
     *                                       simple_index) then
                  print *,'Error: Trailing indices of array ',xhandle,
     *               ' should all be of type simple_index.'
                  call abort_job() 
               endif
            enddo
         endif
      else if (ynindex .gt. xnindex) then
         nindex = xnindex
         do i = nindex+1, ynindex
            index = array_table(c_index_array1+i-1,yhandle)
            if (index_table(c_index_type,index) .ne. simple_index) then
               print *,'Error: Trailing indices of array ',yhandle,
     *            ' should all be of type simple_index.'
               call abort_job() 
            endif
         enddo
      else if (xnindex .eq. ynindex) then
         nindex = xnindex
      endif

c------------------------------------------------------------------------
c   Look up the info pertaining to the destination block.
c   First, determine the number of indices and the relative index of 
c   the destination block.
c------------------------------------------------------------------------

      if (ytype .eq. static_array) then
         daddr = address_table(yhandle)
         indy  = get_index_from_base(daddr, y, 2)
         do i = 1, nindex
            index = array_table(c_index_array1+i-1,yhandle)

c-------------------------------------------------------------------------
c   Get the data ranges to use for the current segment of this index.
c-------------------------------------------------------------------------

            seg   = index_table(c_current_seg, index)
            call get_index_segment(index, seg, segment_table,
     *                        nsegment_table, index_table,
     *                        nindex_table, yrange1(i),
     *                        yrange2(i))

c-------------------------------------------------------------------------
c   Static arrays: use dimensions of entire array.
c-------------------------------------------------------------------------

            ydim1(i) = array_table(c_index_range1+i-1,yhandle)
            ydim2(i) = array_table(c_index_range2+i-1,yhandle)
         enddo
      else

c--------------------------------------------------------------------------
c   Array is not static, managed by blkmgr.
c--------------------------------------------------------------------------

         stack = array_table(c_array_stack,yhandle)
         indy = get_block_data_index(yhandle, y_block, stack, yndx, y)
 
         call get_block_segments(yndx, blkseg)
         do i = 1, nindex

c-------------------------------------------------------------------------
c   Determine the appropriate index to use.
c-------------------------------------------------------------------------

            index = array_table(c_index_array1+i-1,yhandle)
            if (ysub_lookup(i)) then
               superindex = index_table(c_subindex_ptr, index)
               call get_index_segment(superindex, blkseg(i), 
     *                 segment_table, nsegment_table,
     *                 index_table, nindex_table,
     *                 ydim1(i), ydim2(i))

c-------------------------------------------------------------------------
c   Look up the ranges by using the subindex and its current_seg.
c-------------------------------------------------------------------------

               seg = index_table(c_current_seg, index)     
               call get_index_segment(index, seg,
     *                 segment_table, nsegment_table,
     *                 index_table, nindex_table,
     *                 yrange1(i), yrange2(i)) 
            else

c-------------------------------------------------------------------------
c   Get the data ranges to use for the current segment of this index.
c-------------------------------------------------------------------------

               call get_index_segment(index, blkseg(i), segment_table,
     *                        nsegment_table, index_table,
     *                        nindex_table, yrange1(i),
     *                        yrange2(i))

               ydim1(i) = yrange1(i)
               ydim2(i) = yrange2(i) 
            endif
         enddo
      endif
   
c------------------------------------------------------------------------
c   Look up the info pertaining to the source block.
c------------------------------------------------------------------------

      if (xtype .eq. scalar_value) then
         xscalar = .true.
         i = array_table(c_scalar_index, xhandle)
         xval = scalar_table(i)
         indx = 1   ! dummy reference.
      else if (xtype .eq. static_array) then

c--------------------------------------------------------------------------
c   Source block is a static array.
c--------------------------------------------------------------------------

         daddr = address_table(xhandle)
         indx  = get_index_from_base(daddr, x, 2)
         do i = 1, nindex
            index = array_table(c_index_array1+i-1,xhandle)

c-------------------------------------------------------------------------
c   Get the data ranges to use for the current segment of this index.
c-------------------------------------------------------------------------

            seg   = index_table(c_current_seg, index)
            call get_index_segment(index, seg, segment_table,
     *                        nsegment_table, index_table,
     *                        nindex_table, xrange1(i),
     *                        xrange2(i))

c-------------------------------------------------------------------------
c   Static arrays: use dimensions of entire array.
c-------------------------------------------------------------------------

            xdim1(i) = array_table(c_index_range1+i-1,xhandle)
            xdim2(i) = array_table(c_index_range2+i-1,xhandle)
         enddo
      else

c--------------------------------------------------------------------------
c   Array is not static, managed by blkmgr.
c--------------------------------------------------------------------------

         stack = array_table(c_array_stack,xhandle)
         indx = get_block_data_index(xhandle, x_block, stack, xndx, x)
 
         call get_block_segments(xndx, blkseg)
         do i = 1, nindex

c-------------------------------------------------------------------------
c   Determine the appropriate index to use.
c-------------------------------------------------------------------------

            index = array_table(c_index_array1+i-1,xhandle)
            if (xsub_lookup(i)) then
               superindex = index_table(c_subindex_ptr, index)
               call get_index_segment(superindex, blkseg(i), 
     *                 segment_table, nsegment_table,
     *                 index_table, nindex_table,
     *                 xdim1(i), xdim2(i))

c-------------------------------------------------------------------------
c   Look up the ranges by using the subindex and its current_seg.
c-------------------------------------------------------------------------

               seg = index_table(c_current_seg, index)     
               call get_index_segment(index, seg,
     *                 segment_table, nsegment_table,
     *                 index_table, nindex_table,
     *                 xrange1(i), xrange2(i)) 
            else

c-------------------------------------------------------------------------
c   Get the data ranges to use for the current segment of this index.
c-------------------------------------------------------------------------

               call get_index_segment(index, blkseg(i), segment_table,
     *                        nsegment_table, index_table,
     *                        nindex_table, xrange1(i),
     *                        xrange2(i))

               xdim1(i) = xrange1(i)
               xdim2(i) = xrange2(i) 
            endif
         enddo
      endif
   
c      print *,'SLICE_INSERT_BLOCK line ',current_line
c      print *,'Source block: ',xhandle, x_block, xndx
c      print *,'   dim: ',(xdim1(i),xdim2(i),i=1,nindex)
c      print *,'   range: ',(xrange1(i),xrange2(i),i=1,nindex)
c      print *,'   xtype ',xtype,' xscalar ',xscalar
c      if (xscalar) print *,'xval ',xval
c      print *,'Dest block: ',yhandle, y_block, yndx
c      print *,'   dim: ',(ydim1(i),ydim2(i),i=1,nindex)
c      print *,'   range: ',(yrange1(i),yrange2(i),i=1,nindex)
c      print *,'   ytype ',ytype

c-----------------------------------------------------------------------------
c   Perform actual data movement.
c-----------------------------------------------------------------------------

      if (nindex .eq. 2) then
           call slice_insert2(y(indy), ydim1(1), ydim1(2),
     *                         ydim2(1), ydim2(2),
     *                         yrange1(1), yrange1(2),
     *                         yrange2(1), yrange2(2),
     *                         x(indx), xdim1(1), xdim1(2),
     *                         xdim2(1), xdim2(2),
     *                         xscalar, xval)
      else if (nindex .eq. 4) then
           call slice_insert4(y(indy), ydim1(1), ydim1(2), 
     *                         ydim1(3), ydim1(4),
     *                         ydim2(1), ydim2(2), ydim2(3), ydim2(4),
     *                         yrange1(1), yrange1(2), 
     *                         yrange1(3), yrange1(4),
     *                         yrange2(1), yrange2(2),
     *                         yrange2(3), yrange2(4),
     *                         x(indx), xdim1(1), xdim1(2), 
     *                         xdim1(3), xdim1(4),
     *                         xdim2(1), xdim2(2), xdim2(3), xdim2(4),
     *                         xscalar, xval)
      else if (nindex .eq. 6) then
           call slice_insert6(y(indy), ydim1(1), ydim1(2),
     *                         ydim1(3), ydim1(4),ydim1(5), ydim1(6),
     *                         ydim2(1), ydim2(2), ydim2(3), ydim2(4),
     *                         ydim2(5), ydim2(6),
     *                         yrange1(1), yrange1(2),
     *                         yrange1(3), yrange1(4),
     *                         yrange1(5), yrange1(6),
     *                         yrange2(1), yrange2(2),
     *                         yrange2(3), yrange2(4),
     *                         yrange2(5), yrange2(6),
     *                         x(indx), xdim1(1), xdim1(2),
     *                         xdim1(3), xdim1(4), xdim1(5), xdim1(6),
     *                         xdim2(1), xdim2(2), xdim2(3), xdim2(4),
     *                         xdim2(5), xdim2(6),
     *                         xscalar, xval)
      else 
         print *,'Error: No slice/insert available for arrays of ',
     *           'dimension ',nindex
         call abort_job()
      endif

      return
      end

      subroutine slice_insert2(y, ya1, yb1,  
     *                         ya2, yb2, 
     *                         a1, b1, 
     *                         a2, b2, 
     *                         x, xa1, xb1, 
     *                         xa2, xb2, 
     *                         xscalar, xcons)
c----------------------------------------------------------------------------
c   Inserts a slice from x(xa1:xa2,xb1:xb2) into 
c   y(ya1:ya2,yb1:yb2).  If xscalar
c   is true, we insert xcons into the subblock (a1:a2,b1:b2) of 
c   y. 
c
c   The orbital ranges of the dimensions of the sub-block that we are 
c   slicing from/inserting into ar (a1:a2,b1:b2).
c----------------------------------------------------------------------------
      implicit none

      integer ya1, yb1, a1, b1
      integer ya2, yb2, a2, b2
      integer xa1, xb1
      integer xa2, xb2
      double precision y(ya1:ya2,yb1:yb2)
      double precision x(xa1:xa2,xb1:xb2)
      double precision xcons
      logical xscalar

      integer a, b

      if (a1 .lt. ya1 .or. a2 .gt. ya2 .or.
     *    b1 .lt. yb1 .or. b2 .gt. yb2) then
         print *,'SLICE_INSERT_BLOCK: Out of range error.'
         print *,'Destination block indices: ',ya1,ya2,yb1,yb2
         print *,'Sub-block target indices: ',a1,a2,b1,b2
         call abort_job()
      endif

      if (xscalar) then

c----------------------------------------------------------------------------
c   Scalar insertion.
c----------------------------------------------------------------------------

         do b = b1, b2
         do a = a1, a2
            y(a,b) = xcons
         enddo
         enddo
      else

c----------------------------------------------------------------------------
c   Normal (possible) insertion from a (possible) slice.
c----------------------------------------------------------------------------

         if (a1 .lt. xa1 .or. a2 .gt. xa2 .or.
     *       b1 .lt. xb1 .or. b2 .gt. xb2) then
            print *,'SLICE_INSERT_BLOCK: Out of range error.'
            print *,'Source block indices: ',xa1,xa2,xb1,xb2
            print *,'Sub-block target indices: ',a1,a2,b1,b2
            call abort_job()
         endif

         do b = b1, b2  
         do a = a1, a2  
            y(a,b) = x(a,b)
         enddo
         enddo
      endif
 
      return
      end

      subroutine slice_insert4(y, ya1, yb1, yc1, yd1, 
     *                         ya2, yb2, yc2, yd2,
     *                         a1, b1, c1, d1,
     *                         a2, b2, c2, d2,
     *                         x, xa1, xb1, xc1, xd1, 
     *                         xa2, xb2, xc2, xd2,
     *                         xscalar, xcons)
c----------------------------------------------------------------------------
c   Inserts a slice from x(xa1:xa2,xb1:xb2,xc1:xc2,xd1:xd2) into 
c   y(ya1:ya2,yb1:yb2,yc1:yc2,yd1:yd2).  If xscalar
c   is true, we insert xcons into the subblock (a1:a2,b1:b2,c1:c2,d1:d2) of 
c   y. 
c
c   The orbital ranges of the dimensions of the sub-block that we are 
c   slicing from/inserting into ar (a1:a2,b1:b2,c1:c2,d1:d2).
c----------------------------------------------------------------------------
      implicit none

      integer ya1, yb1, yc1, yd1, a1, b1, c1, d1
      integer ya2, yb2, yc2, yd2, a2, b2, c2, d2
      integer xa1, xb1, xc1, xd1
      integer xa2, xb2, xc2, xd2
      double precision y(ya1:ya2,yb1:yb2,yc1:yc2,yd1:yd2)
      double precision x(xa1:xa2,xb1:xb2,xc1:xc2,xd1:xd2)
      double precision xcons
      logical xscalar

      integer a, b, c, d

      if (a1 .lt. ya1 .or. a2 .gt. ya2 .or.
     *    b1 .lt. yb1 .or. b2 .gt. yb2 .or.
     *    c1 .lt. yc1 .or. c2 .gt. yc2 .or.
     *    d1 .lt. yd1 .or. d2 .gt. yd2) then
         print *,'SLICE_INSERT_BLOCK: Out of range error.'
         print *,'Destination block indices: ',ya1,ya2,yb1,yb2,
     *            yc1, yc2, yd1, yd2
         print *,'Sub-block target indices: ',a1,a2,b1,b2,c1,c2,d1,d2
         call abort_job()
      endif

      if (xscalar) then

c----------------------------------------------------------------------------
c   Scalar insertion.
c----------------------------------------------------------------------------

         do d = d1, d2
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
            y(a,b,c,d) = xcons
         enddo
         enddo
         enddo
         enddo
      else

c----------------------------------------------------------------------------
c   Normal (possible) insertion from a (possible) slice.
c----------------------------------------------------------------------------

         if (a1 .lt. xa1 .or. a2 .gt. xa2 .or.
     *       b1 .lt. xb1 .or. b2 .gt. xb2 .or.
     *       c1 .lt. xc1 .or. c2 .gt. xc2 .or.
     *       d1 .lt. xd1 .or. d2 .gt. xd2) then
            print *,'SLICE_INSERT_BLOCK: Out of range error.'
            print *,'Source block indices: ',xa1,xa2,xb1,xb2,
     *            xc1, xc2, xd1, xd2
            print *,'Sub-block target indices: ',a1,a2,b1,b2,c1,c2,d1,d2
            call abort_job()
         endif

         do d = d1, d2  
         do c = c1, c2  
         do b = b1, b2  
         do a = a1, a2  
            y(a,b,c,d) = x(a,b,c,d)
         enddo
         enddo
         enddo
         enddo
      endif
 
      return
      end

      subroutine slice_insert6(y, ya1, yb1, yc1, yd1,ye1,yf1, 
     *                         ya2, yb2, yc2, yd2, ye2, yf2,
     *                         a1, b1, c1, d1, e1, f1,
     *                         a2, b2, c2, d2, e2,f2,
     *                         x, xa1, xb1, xc1, xd1, xe1, xf1,
     *                         xa2, xb2, xc2, xd2, xe2, xf2,
     *                         xscalar, xcons)
c----------------------------------------------------------------------------
c   Inserts a slice from x(xa1:xa2,xb1:xb2,xc1:xc2,xd1:xd2,xe1:xe2,xf1:xf2) into 
c   y(ya1:ya2,yb1:yb2,yc1:yc2,yd1:yd2,ye1:ye2,yf1:yf2).  If xscalar
c   is true, we insert xcons into the subblock 
c   (a1:a2,b1:b2,c1:c2,d1:d2,e1:e2,f1:f2) of y.
c
c   The orbital ranges of the dimensions of the sub-block that we are 
c   slicing from/inserting into ar (a1:a2,b1:b2,c1:c2,d1:d2,e1:e2,f1:f2).
c----------------------------------------------------------------------------
      implicit none

      integer ya1, yb1, yc1, yd1, ye1,yf1, a1, b1, c1, d1, e1, f1
      integer ya2, yb2, yc2, yd2, ye2, yf2, a2, b2, c2, d2, e2,f2
      integer xa1, xb1, xc1, xd1, xe1, xf1
      integer xa2, xb2, xc2, xd2, xe2, xf2 
      double precision 
     *   y(ya1:ya2,yb1:yb2,yc1:yc2,yd1:yd2,ye1:ye2,yf1:yf2), 
     *   x(xa1:xa2,xb1:xb2,xc1:xc2,xd1:xd2,xe1:xe2,xf1:xf2)
      double precision xcons
      logical xscalar

      integer a, b, c, d, e, f

      if (a1 .lt. ya1 .or. a2 .gt. ya2 .or.
     *    b1 .lt. yb1 .or. b2 .gt. yb2 .or.
     *    c1 .lt. yc1 .or. c2 .gt. yc2 .or.
     *    d1 .lt. yd1 .or. d2 .gt. yd2 .or.
     *    e1 .lt. ye1 .or. e2 .gt. ye2 .or.
     *    f1 .lt. yf1 .or. f2 .gt. yf2) then
         print *,'SLICE_INSERT_BLOCK: Out of range error.'
         print *,'Destination block indices: ',ya1,ya2,yb1,yb2,
     *            yc1, yc2, yd1, yd2, ye1, ye2, yf1, yf2
         print *,'Sub-block target indices: ',a1,a2,b1,b2,c1,c2,d1,d2,
     *            e1, e2, f1, f2
         call abort_job()
      endif

      if (xscalar) then

c----------------------------------------------------------------------------
c   Scalar insertion.
c----------------------------------------------------------------------------

         do f = f1, f2
         do e = e1, e2
         do d = d1, d2
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
            y(a,b,c,d,e,f) = xcons
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      else

c----------------------------------------------------------------------------
c   Normal (possible) insertion from a (possible) slice.
c----------------------------------------------------------------------------

         if (a1 .lt. xa1 .or. a2 .gt. xa2 .or.
     *       b1 .lt. xb1 .or. b2 .gt. xb2 .or.
     *       c1 .lt. xc1 .or. c2 .gt. xc2 .or.
     *       d1 .lt. xd1 .or. d2 .gt. xd2 .or.
     *       e1 .lt. xe1 .or. e2 .gt. xe2 .or.
     *       f1 .lt. xf1 .or. f2 .gt. xf2) then
            print *,'SLICE_INSERT_BLOCK: Out of range error.'
            print *,'Source block indices: ',xa1,xa2,xb1,xb2,
     *            xc1, xc2, xd1, xd2, xe1, xe2, xf1, xf2
            print *,'Sub-block target indices: ',
     *          a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2
            call abort_job()
         endif

         do f = f1, f2
         do e = e1, e2
         do d = d1, d2  
         do c = c1, c2  
         do b = b1, b2  
         do a = a1, a2  
            y(a,b,c,d,e,f) = x(a,b,c,d,e,f)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      endif
 
      return
      end

