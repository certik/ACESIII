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
      subroutine sum_blocks(x1handle, x1_block, x1_ndx, 
     *                      x2handle, x2_block, x2_ndx,
     *                      yhandle, y_block, y_ndx, opcode,
     *                      array_table, narray_table, 
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table, 
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
c--------------------------------------------------------------------------
c   Forms the sum of blocks x1 and x2, resulting in block y.
c   x1, x2, and y are the array_table entries for their respective arrays.
c   x1_block, x2_block, and y_block are the actual block indices of the
c   blocks to be summed.
c
c   If the opcode is a subtract_op, this routine forms the difference
c   of x1-x2, storing the result in y.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer nsegment_table 
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      integer x1handle, x2handle, yhandle
      integer x1_block, x2_block, y_block 
      integer x1_ndx, x2_ndx, y_ndx
      integer opcode, flopcount

      integer n1, n2, ny
      integer nw, iseg
      integer*8 indy, indx1, indx2, daddr
      integer*8 get_block_data_index
      integer stack
      integer*8 get_index_from_base
      integer i, j, k
      integer ind(mx_array_index,3)
      integer range1(mx_array_index), range2(mx_array_index)
      integer dim1(mx_array_index,3), dim2(mx_array_index,3)
      integer blkmgr_blocksize
      integer x1type, x2type, ytype

      double precision y(1), x1(1), x2(2)
#ifdef ALTIX
      pointer (ptry, y)
      pointer (ptrx1, x1)
      pointer (ptrx2, x2)
#else
      common y, x1, x2
#endif

      x1type = array_table(c_array_type, x1handle)
      x2type = array_table(c_array_type, x2handle)
      ytype = array_table(c_array_type, yhandle)
      
#ifdef ALTIX
      ptry = dshptr
      ptrx1 = dshptr
      ptrx2 = dshptr
#endif

c--------------------------------------------------------------------------
c   Verify no operands are scalars.
c--------------------------------------------------------------------------

      if (x1type .eq. scalar_value .or.
     *    x2type .eq. scalar_value .or.
     *    ytype .eq. scalar_value) then
         print *,'Error in sum_blocks: One or more scalar operands.'
         print *,'Operands are ',x1handle,x2handle,yhandle
      endif
 
      n1 = array_table(c_nindex, x1handle)
      n2 = array_table(c_nindex, x2handle)
      ny = array_table(c_nindex, yhandle)

      flopcount = 0

      do j = 1, 3
      do i = 1, mx_array_index
        ind(i,j) = 0
      enddo
      enddo
     
c-----------------------------------------------------------------------
c   Match the operand indices with that of the result.
c-----------------------------------------------------------------------

      do i = 1, ny
         ind(i,1) = array_table(c_index_array1+i-1,yhandle)
         
         do j = 1, n1
            if (ind(i,1) .eq. array_table(c_index_array1+j-1,x1handle)) 
     *            ind(i,2) = array_table(c_index_array1+j-1,x1handle)
         enddo

         do j = 1, n2
            if (ind(i,1) .eq. array_table(c_index_array1+j-1,x2handle))
     *            ind(i,3) = array_table(c_index_array1+j-1,x2handle)
         enddo
      enddo

c--------------------------------------------------------------------------
c   All indices must match for a summation operation.
c--------------------------------------------------------------------------

      do i = 1, ny
         if (ind(i,1) .ne. ind(i,2) .or.
     *       ind(i,1) .ne. ind(i,3)) then
            print *,'Error: Summation of blocks, non-matching indices.'
            print *,'Indices:'
            do j = 1, ny
               print *,'   ',ind(j,1),ind(j,2),ind(j,3)
            enddo
            call abort_job()
         endif
      enddo

c---------------------------------------------------------------------------
c   Compute the size of the block.
c---------------------------------------------------------------------------

      nw = 1
      do i = 1, ny
         iseg = index_table(c_current_seg,ind(i,1))
         call get_index_segment(ind(i,1), iseg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, range1(i), range2(i))
         nw = nw * (range2(i)-range1(i)+1)
      enddo

      if (ytype .eq. static_array) then
         daddr = address_table(yhandle)
         indy  = get_index_from_base(daddr, y, 2)
         do i = 1, ny
            dim1(i,1) = array_table(c_index_range1+i-1,yhandle)
            dim2(i,1) = array_table(c_index_range2+i-1,yhandle)
         enddo
      else
         stack = array_table(c_array_stack,yhandle)
         indy  = get_block_data_index(yhandle, y_block, stack,
     *                                y_ndx, y) 
         do i = 1, ny
            dim1(i,1) = range1(i)
            dim2(i,1) = range2(i)
         enddo
      endif

      if (x1type .eq. static_array) then
         daddr = address_table(x1handle)
         indx1 = get_index_from_base(daddr, x1, 2)
         do i = 1, n1
            dim1(i,2) = array_table(c_index_range1+i-1,x1handle)
            dim2(i,2) = array_table(c_index_range2+i-1,x1handle)
         enddo
      else
         stack = array_table(c_array_stack,x1handle)
         indx1 = get_block_data_index(x1handle, x1_block, stack, x1_ndx,
     *                                      x1) 
         do i = 1, n1
            dim1(i,2) = range1(i)
            dim2(i,2) = range2(i)
         enddo
      endif

      if (x2type .eq. static_array) then
         daddr = address_table(x2handle)
         indx2 = get_index_from_base(daddr, x2, 2)
         do i = 1, n2
            dim1(i,3) = array_table(c_index_range1+i-1,x2handle)
            dim2(i,3) = array_table(c_index_range2+i-1,x2handle)
         enddo
      else
         stack = array_table(c_array_stack,x2handle)
         indx2 = get_block_data_index(x2handle, x2_block, stack,
     *                                x2_ndx, x2) 
         do i = 1, n2
            dim1(i,3) = range1(i)
            dim2(i,3) = range2(i)
         enddo
      endif

c      if (ytype .ne. static_array .and. 
c     *    x1type .ne. static_array .and.
c     *    x2type .ne. static_array) then

c---------------------------------------------------------------------------
c   Sum the blocks.
c---------------------------------------------------------------------------

c         call add_sub_blocks(y(indy), x1(indx1), x2(indx2), nw,
c     *                          opcode)
c      else

c----------------------------------------------------------------------------
c   Special code to handle static array dimensions.
c----------------------------------------------------------------------------

         if (ny .eq. 4) then
            call sum_static_blocks4(y(indy),range1(1),range2(1),
     *             range1(2),range2(2),range1(3),range2(3),
     *             range1(4),range2(4),dim1(1,1),dim2(1,1),
     *             dim1(2,1),dim2(2,1),dim1(3,1),dim2(3,1),
     *             dim1(4,1),dim2(4,1), x1(indx1), 
     *             dim1(1,2),dim2(1,2),dim1(2,2),dim2(2,2),
     *             dim1(3,2),dim2(3,2),dim1(4,2),dim2(4,2),
     *             x2(indx2),
     *             dim1(1,3),dim2(1,3),dim1(2,3),dim2(2,3),
     *             dim1(3,3),dim2(3,3),dim1(4,3),dim2(4,3),
     *             opcode)
         else if (ny .eq. 2) then
            call sum_static_blocks2(y(indy),range1(1),range2(1),
     *             range1(2),range2(2),
     *             dim1(1,1),dim2(1,1),dim1(2,1),dim2(2,1),
     *             x1(indx1),
     *             dim1(1,2),dim2(1,2),dim1(2,2),dim2(2,2),
     *             x2(indx2),
     *             dim1(1,3),dim2(1,3),dim1(2,3),dim2(2,3),
     *             opcode)
         else if (ny .eq. 6) then
            call sum_static_blocks6(y(indy),range1(1),range2(1),
     *             range1(2),range2(2),range1(3),range2(3),
     *             range1(4),range2(4),range1(5),range2(5),
     *             range1(6),range2(6),
     *             dim1(1,1),dim2(1,1),
     *             dim1(2,1),dim2(2,1),dim1(3,1),dim2(3,1),
     *             dim1(4,1),dim2(4,1), 
     *             dim1(5,1),dim2(5,1), dim1(6,1),dim2(6,1), 
     *             x1(indx1),
     *             dim1(1,2),dim2(1,2),dim1(2,2),dim2(2,2),
     *             dim1(3,2),dim2(3,2),dim1(4,2),dim2(4,2),
     *             dim1(5,2),dim2(5,2),dim1(6,2),dim2(6,2),
     *             x2(indx2),
     *             dim1(1,3),dim2(1,3),dim1(2,3),dim2(2,3),
     *             dim1(3,3),dim2(3,3),dim1(4,3),dim2(4,3),
     *             dim1(5,3),dim2(5,3),dim1(6,3),dim2(6,3),
     *             opcode)
         else if (ny .eq. 3) then 
            call sum_static_blocks3(y(indy),range1(1),range2(1),
     *             range1(2),range2(2),range1(3),range2(3),
     *             dim1(1,1),dim2(1,1),
     *             dim1(2,1),dim2(2,1),dim1(3,1),dim2(3,1),
     *             x1(indx1),
     *             dim1(1,2),dim2(1,2),dim1(2,2),dim2(2,2),
     *             dim1(3,2),dim2(3,2),
     *             x2(indx2),
     *             dim1(1,3),dim2(1,3),dim1(2,3),dim2(2,3),
     *             dim1(3,3),dim2(3,3),
     *             opcode)
         else
            print *,'Error: Cannot handle static arrays of dimension',
     *           ' other than 2,3, 4, or 6.'
            print *,'Arrays have dimension of ',ny
            call abort_job()
         endif
c      endif

      flopcount = nw

      return
      end

      subroutine sum_static_blocks6(y, a1,a2,b1,b2,c1,c2,d1,d2,
     *                              e1,e2,f1,f2,
     *                              na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                              ny51,ny52,ny61,ny62,
     *                              x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *                              nx1_51,nx1_52,nx1_61,nx1_62,
     *                              x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                              nx2_51,nx2_52,nx2_61,nx2_62,
     *                              opcode)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2
      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ny51,ny52,ny61,ny62
      integer ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,nx1_51,nx1_52,nx1_61,
     *        nx1_62
      integer ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,nx2_51,nx2_52,nx2_61,
     *        nx2_62
      double precision 
     *    y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,ny51:ny52,ny61:ny62)
      double precision 
     *  x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2,nx1_51:nx1_52,nx1_61:nx1_62)
      double precision 
     *  x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2,nx2_51:nx2_52,nx2_61:nx2_62)
      integer opcode, a,b,c,d,e,f

      if (opcode .eq. subtract_op) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,b,c,d,e,f) - x2(a,b,c,d,e,f)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      else if (opcode .eq. sum_op) then
         do f = f1,f2
         do e = e1,e2
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d,e,f) = x1(a,b,c,d,e,f) + x2(a,b,c,d,e,f) 
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         print *,'Error: Invalid opcode in sum_static_blocks6: ',
     *      opcode
         call abort_job()
      endif

      return
      end

      subroutine sum_static_blocks4(y, a1,a2,b1,b2,c1,c2,d1,d2,
     *                              na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                              x1,ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2,
     *                              x2,ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2,
     *                              opcode)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2
      integer ne1,ne2,nf1,nf2,ng1,ng2,nh1,nh2
      integer ni1,ni2,nj1,nj2,nk1,nk2,nl1,nl2
      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2,nh1:nh2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2,nl1:nl2) 
      integer opcode, a,b,c,d

      if (opcode .eq. subtract_op) then
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d) = x1(a,b,c,d) - x2(a,b,c,d)
         enddo
         enddo
         enddo
         enddo
      else if (opcode .eq. sum_op) then
         do d = d1,d2
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c,d) = x1(a,b,c,d) + x2(a,b,c,d)
         enddo
         enddo
         enddo
         enddo
      else
         print *,'Error: Invalid opcode in sum_static_blocks4: ',
     *      opcode
         call abort_job()
      endif

      return
      end

      subroutine sum_static_blocks3(y, a1,a2,b1,b2,c1,c2,
     *                              na1,na2,nb1,nb2,nc1,nc2,
     *                              x1,ne1,ne2,nf1,nf2,ng1,ng2,
     *                              x2,ni1,ni2,nj1,nj2,nk1,nk2,
     *                              opcode)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,c1,c2
      integer na1,na2,nb1,nb2,nc1,nc2
      integer ne1,ne2,nf1,nf2,ng1,ng2
      integer ni1,ni2,nj1,nj2,nk1,nk2
      double precision y(na1:na2,nb1:nb2,nc1:nc2)
      double precision x1(ne1:ne2,nf1:nf2,ng1:ng2)
      double precision x2(ni1:ni2,nj1:nj2,nk1:nk2) 
      integer opcode, a,b,c,d

      if (opcode .eq. subtract_op) then
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c) = x1(a,b,c) - x2(a,b,c)
         enddo
         enddo
         enddo
      else if (opcode .eq. sum_op) then
         do c = c1,c2
         do b = b1,b2
         do a = a1,a2
            y(a,b,c) = x1(a,b,c) + x2(a,b,c)
         enddo
         enddo
         enddo
      else
         print *,'Error: Invalid opcode in sum_static_blocks3: ',
     *      opcode
         call abort_job()
      endif

      return
      end

      subroutine sum_static_blocks2(y, a1,a2,b1,b2,
     *                              na1,na2,nb1,nb2,
     *                              x1,ne1,ne2,nf1,nf2,
     *                              x2,ni1,ni2,nj1,nj2,
     *                              opcode)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2
      integer na1,na2,nb1,nb2
      integer ne1,ne2,nf1,nf2
      integer ni1,ni2,nj1,nj2
      double precision y(na1:na2,nb1:nb2)
      double precision x1(ne1:ne2,nf1:nf2)
      double precision x2(ni1:ni2,nj1:nj2)
      integer opcode, a,b

      if (opcode .eq. subtract_op) then
         do b = b1,b2
         do a = a1,a2
            y(a,b) = x1(a,b) - x2(a,b)
         enddo
         enddo
      else if (opcode .eq. sum_op) then
         do b = b1,b2
         do a = a1,a2
            y(a,b) = x1(a,b) + x2(a,b)
         enddo
         enddo
      else
         print *,'Error: Invalid opcode in sum_static_blocks2: ',
     *      opcode
         call abort_job()
      endif

      return
      end

      subroutine add_sub_blocks(y, x1, x2, n, opcode)
      implicit none
      include 'interpreter.h'

      integer n, opcode
      double precision y(n), x1(n), x2(n)

      integer i

      if (opcode .eq. sum_op) then
         do i = 1, n
            y(i) = x1(i) + x2(i)
         enddo    
      else if (opcode .eq. subtract_op) then
         do i = 1, n
            y(i) = x1(i) - x2(i)
         enddo    
      else
         print *,'Error in sum_blocks: Invalid opcode ',opcode
         call abort_job()
      endif

      return
      end
