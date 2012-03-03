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
c   is performed.
      subroutine contract_blocks(x1handle, x1_block, x1_ndx,
     *                      x2handle, x2_block, x2_ndx,
     *                      yhandle, y_block, y_ndx, cind, ncind,
     *                      array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, block_map_table,
     *                      flopcount)
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'machine_types.h'
      include 'saved_data.h'
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
      integer block_map_table(*)
      integer ncind, cind(ncind)
      integer x1handle, x2handle, yhandle
      integer x1_block, x2_block, y_block 
      integer x1_ndx, x2_ndx, y_ndx
      integer flopcount
      integer iseq

      integer n1, n2, ny
      integer nw, iseg, range1, range2
      integer*8 iy, ix1, ix2, daddr
      integer*8 get_block_index, get_index_from_base
      integer ystack, x1stack, x2stack
      integer i, j, k
      integer ind(4,3)
      integer x1_type, x2_type, ytype
      integer nwblock
      double precision value

      double precision y(1), x1(1), x2(1)
#ifdef ALTIX
      pointer (yptr, y)
      pointer (x1ptr, x1)
      pointer (x2ptr, x2)
      pointer (scr_ptr, scr)
#else
      common y, x1, x2
#endif

      integer indy(mx_array_index), segy(mx_array_index), 
     *        ya(mx_array_index), yb(mx_array_index)
      integer indx1(mx_array_index), segx1(mx_array_index), 
     *        x1a(mx_array_index), x1b(mx_array_index)
      integer indx2(mx_array_index), segx2(mx_array_index), 
     *        x2a(mx_array_index), x2b(mx_array_index)
      integer nya(mx_array_index), nyb(mx_array_index), 
     *        nx1a(mx_array_index), nx1b(mx_array_index), 
     *        nx2a(mx_array_index), nx2b(mx_array_index)

      integer handle1, handle2, handle3, istat
      integer allocate_scratch_block
      integer*8 iscr1, iscr2, iscr3
      integer*8 c_loc64, n
      double precision scr(1)
      save handle1, handle2, handle3
      save iscr1, iscr2, iscr3

#ifdef ALTIX
      yptr = dshptr
      x1ptr = dshptr
      x2ptr = dshptr
      scr_ptr = dshptr
#endif

      if (first_contract_blocks) then
         first_contract_blocks = .false.
 
c------------------------------------------------------------------------
c   Allocate 3 scratch blocks for the 2-index contraction routines.
c------------------------------------------------------------------------

         istat = allocate_scratch_block(scr, iscr1, handle1,
     *                           array_table, narray_table,
     *                           index_table, nindex_table, 
     *                           block_map_table)
         if (istat .ne. 0) then
            print *,'Error: Unable to allocate scratch blocks.'
            call abort_job()
         endif

         istat = allocate_scratch_block(scr, iscr2, handle2,
     *                           array_table, narray_table,
     *                           index_table, nindex_table,
     *                           block_map_table)
         if (istat .ne. 0) then
            print *,'Error: Unable to allocate scratch blocks.'
            call abort_job()
         endif

         istat = allocate_scratch_block(scr, iscr3, handle3,
     *                           array_table, narray_table,
     *                           index_table, nindex_table,
     *                           block_map_table)
         if (istat .ne. 0) then
            print *,'Error: Unable to allocate scratch blocks.'
            call abort_job()
         endif
      endif

c      print *,'CONTRACT_BLOCKS: x1 = ',x1handle,' x1_block = ',x1_block,
c     *        ' x2 = ',x2handle,' x2_block = ',x2_block,' y = ',yhandle,
c     *        ' y_block = ',y_block,' indices = ',(cind(i),i=1,ncind)

      if ((ytype.eq.scalar_value) .and. (x1_type.eq.scalar_value) .and. 
     *    (x2_type.eq.scalar_value)) then 
        write(6,*) 'Scalar only mult being skipped ' 
        return 
      endif 

      flopcount = 0

c------------------------------------------------------------------------
c   Chcek for any scalar operands.
c------------------------------------------------------------------------

      x1_type = array_table(c_array_type,x1handle)
      x2_type = array_table(c_array_type,x2handle)
      ytype   = array_table(c_array_type,yhandle)

      x1stack = array_table(c_array_stack,x1handle)
      x2stack = array_table(c_array_stack,x2handle)
      ystack  = array_table(c_array_stack,yhandle)

      if (x1_type .eq. scalar_value) then
         value = scalar_table(array_table(c_scalar_index,x1handle))

         if (x2_type .eq. static_array .or.
     *       ytype   .eq. static_array) then
            call multiply_static_block(x2handle, x2_block, x2_ndx,
     *                value, yhandle, y_block, y_ndx,
     *                array_table, narray_table,
     *                index_table, nindex_table, 
     *                segment_table, nsegment_table,
     *                address_table)
         else
            call get_actual_blocksize(x2handle, x2_block, x2_ndx,
     *              array_table, narray_table,
     *              index_table, nindex_table,
     *              segment_table, nsegment_table, nwblock)
            
            call multiply_block(x2handle, x2_block, x2stack, x2_ndx,  
     *                       value, yhandle, y_block, ystack, 
     *                       y_ndx, nwblock)
         endif
         return
      else if (x2_type .eq. scalar_value) then
         value = scalar_table(array_table(c_scalar_index,x2handle))
         if (x1_type .eq. static_array .or.
     *       ytype   .eq. static_array) then
            call multiply_static_block(x1handle, x1_block, x1_ndx,
     *                value, yhandle, y_block, y_ndx,
     *                array_table, narray_table,
     *                index_table, nindex_table,
     *                segment_table, nsegment_table,
     *                address_table)
         else
            call get_actual_blocksize(x1handle, x1_block, x1_ndx,
     *              array_table, narray_table,
     *              index_table, nindex_table,
     *              segment_table, nsegment_table, nwblock)
            call multiply_block(x1handle, x1_block, x1stack, x1_ndx,  
     *                       value, yhandle, y_block, ystack,
     *                       y_ndx, nwblock)
         endif
         return
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

      if (ytype .eq. scalar_value) go to 100

      do i = 1, ny
         indy(i) = array_table(c_index_array1+i-1, yhandle)
         segy(i) = index_table(c_current_seg,indy(i))
         call get_index_segment(indy(i), segy(i), segment_table,
     *                     nsegment_table, index_table,
     *                     nindex_table, ya(i), yb(i))  
      enddo

      if (ytype .eq. static_array) then
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
         iy = get_block_index(yhandle, y_block, ystack, 
     *                        y_ndx, y, .true.)
         do i = 1, ny
            nya(i) = ya(i)
            nyb(i) = yb(i)
         enddo
      endif

  100 continue

      do i = 1, n1
         indx1(i) = array_table(c_index_array1+i-1, x1handle)
         segx1(i) = index_table(c_current_seg,indx1(i))
         call get_index_segment(indx1(i), segx1(i), segment_table,
     *                     nsegment_table, index_table,
     *                     nindex_table, x1a(i), x1b(i))  
      enddo

      if (x1_type .eq. static_array) then
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
         ix1 = get_block_index(x1handle, x1_block, x1stack,
     *                         x1_ndx, x1, .true.)
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
      enddo

      if (x2_type .eq. static_array) then
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
         ix2 = get_block_index(x2handle, x2_block, x2stack,
     *                         x2_ndx, x2, .true.)
         do i = 1, n2
            nx2a(i) = x2a(i)
            nx2b(i) = x2b(i)
         enddo
      endif

c      print *,'Contract_blocks: '
c      print *,'indy: ',(indy(j),j=1,ny)
c      print *,'indx1: ',(indx1(j),j=1,n1)
c      print *,'indx2: ',(indx2(j),j=1,n2)
c      print *,'cind: ',(cind(j),j=1,ncind)

c----------------------------------------------------------------------------
c   Determine the array in which the first index of "y" is found.
c   This determines the order of the operand arrays in some of the
c   calls to low-order contraction routines.
c----------------------------------------------------------------------------

      iseq = 2
      do i = 1, n1
         if (indy(1) .eq. indx1(i)) iseq = 1
      enddo

      if (ny .eq. 2) then
         if (n1 .eq. 4 .and. n2 .eq. 2) then
            if (trace .and.
     *          and(tracelevel,contraction_trace) .ne. 0) 
     *          print *,'cwork2422 at instruction ',current_line

            call cwork2422(y(iy), nya(1), nyb(1), nya(2), nyb(2),
     *                     ya(1), yb(1), ya(2), yb(2),
     *                     indy(1), indy(2), x1(ix1), 
     *                     nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *                     nx1a(3), nx1b(3), nx1a(4), nx1b(4),
     *                     x1a(1),x1b(1),
     *                     x1a(2),x1b(2), x1a(3),x1b(3),x1a(4),x1b(4),
     *                     indx1,
     *                     x2(ix2),nx2a(1), nx2b(1), nx2a(2), nx2b(2), 
     *                     x2a(1), x2b(1), x2a(2), x2b(2), 
     *                     indx2(1), indx2(2),   
     *                     scr(iscr1), scr(iscr2), flopcount)
         else if (n1 .eq. 2 .and. n2 .eq. 4) then

c---------------------------------------------------------------------------
c   Reverse the args, call cwork_2422.
c---------------------------------------------------------------------------

            call cwork2422(y(iy), nya(1), nyb(1), nya(2), nyb(2),
     *                     ya(1), yb(1), ya(2), yb(2),
     *                     indy(1), indy(2), x2(ix2),
     *                     nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *                     nx2a(3), nx2b(3), nx2a(4), nx2b(4),
     *                     x2a(1),x2b(1),
     *                     x2a(2),x2b(2), x2a(3),x2b(3),x2a(4),x2b(4),
     *                     indx2,
     *                     x1(ix1),nx1a(1), nx1b(1), nx1a(2), nx1b(2),
     *                     x1a(1), x1b(1), x1a(2), x1b(2),
     *                     indx1(1), indx1(2),
     *                     scr(iscr1), scr(iscr2), flopcount)
         else if (n1 .eq. 2 .and. n2 .eq. 2 .and. ncind .eq. 1) then
            if (trace .and.
     *          and(tracelevel, contraction_trace) .ne. 0) 
     *          print *,'cwork2221 at instruction ',current_line,
     *          ' iseq ',iseq
            if (iseq .eq. 1) then
               call cwork2221(y(iy), nya(1),nyb(1), nya(2),nyb(2),
     *              ya(1),yb(1), ya(2),yb(2),
     *              indy(1),indy(2),
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              x1a(1),x1b(1),x1a(2),x1b(2),
     *              indx1, x2(ix2), nx2a(1),
     *              nx2b(1),nx2a(2),nx2b(2),
     *              x2a(1), x2b(1), x2a(2),x2b(2), indx2(1),indx2(2), 
     *              scr(iscr1),scr(iscr2))
            else
               call cwork2221(y(iy), nya(1),nyb(1), nya(2),nyb(2),
     *              ya(1),yb(1), ya(2),yb(2),
     *              indy(1),indy(2),
     *              x2(ix2), nx2a(1),
     *              nx2b(1),nx2a(2),nx2b(2),
     *              x2a(1), x2b(1), x2a(2),x2b(2), indx2, 
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              x1a(1),x1b(1),x1a(2),x1b(2),indx1(1),indx1(2),
     *              scr(iscr1),scr(iscr2))
            endif

         else if (n1 .eq. 4 .and. n2 .eq. 4 .and. ncind .eq. 3) then
            if (trace .and.
     *          and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork2443 at instruction ',current_line
            call cwork2443(y(iy),nya(1),nyb(1),nya(2),nyb(2),
     *                     ya(1),yb(1),ya(2),yb(2),indy(1), indy(2),
     *                     x1(ix1),nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *                     nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *                     x1a(1),x1b(1),x1a(2),x1b(2),
     *                     x1a(3),x1b(3),x1a(4),x1b(4),indx1,
     *                     x2(ix2),nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *                     nx2a(3),nx2b(3),nx2a(4),nx2b(4),
     *                     x2a(1), x2b(1), x2a(2),x2b(2),
     *                     x2a(3), x2b(3), x2a(4),x2b(4),indx2,
     *                     flopcount,scr(iscr1),scr(iscr2))

         else if (n1 .eq. 6 .and. n2 .eq. 4 .and. ncind .eq. 4) then
c           return 
            call cwork2644(y(iy),nya(1),nyb(1),nya(2),nyb(2),
     *                     ya(1),yb(1),ya(2),yb(2),indy(1), indy(2),

     *                     x1(ix1),nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *                     nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *                     nx1a(5),nx1b(5),nx1a(6),nx1b(6),
     *                     x1a(1),x1b(1),x1a(2),x1b(2),
     *                     x1a(3),x1b(3),x1a(4),x1b(4),
     *                     x1a(5),x1b(5),x1a(6),x1b(6),indx1,

     *                     x2(ix2),nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *                     nx2a(3),nx2b(3),nx2a(4),nx2b(4),
     *                     x2a(1), x2b(1), x2a(2),x2b(2),
     *                     x2a(3), x2b(3), x2a(4),x2b(4),indx2,
     *                     cind, 
     *                     flopcount,scr(iscr1),scr(iscr2),scr(iscr3))


         else
            print *,'Error: Invalid contraction'
            print *,'ny, n1, n2, ncind = ',ny, n1, n2, ncind
            call abort_job()
         endif
      else if (ny .eq. 1) then
         if ((n1 .eq. 3) .and. (n2 .eq. 2)) then

            call cwork1322(y(iy), nya(1), nyb(1), ya(1), yb(1),
     *                     indy(1), x1(ix1),
     *                     nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *                     nx1a(3), nx1b(3),
     *                     x1a(1),x1b(1),
     *                     x1a(2),x1b(2), x1a(3),x1b(3),
     *                     indx1,
     *                     x2(ix2),nx2a(1), nx2b(1), nx2a(2), nx2b(2),
     *                     x2a(1), x2b(1), x2a(2), x2b(2),
     *                     indx2(1), indx2(2),
     *                     scr(iscr1), scr(iscr2), flopcount)
         endif
      else if (ny .eq. 3) then
         if ((n1 .eq. 3) .and. (n2 .eq. 2) .and. (ncind .eq. 1)) then
            call cwork3321(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3),  
     *              ya(1),yb(1),ya(2),yb(2),ya(3),yb(3),
     *              indy(1),indy(2),indy(3),
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              indx1, x2(ix2), nx2a(1), nx2b(1),
     *              nx2a(2), nx2b(2), x2a(1),
     *              x2b(1),x2a(2),x2b(2), indx2(1),indx2(2), cind, 
     *              flopcount, scr(iscr1), scr(iscr2))
         endif

      else if (ny .eq. 4) then

         if ((n1 .eq. 3) .and. (n2 .eq. 3) .and. (ncind .eq. 1)) then
            call cwork4331(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3),nya(4),nyb(4),  
     *              ya(1),yb(1),ya(2),yb(2),ya(3),yb(3),ya(4),yb(4),
     *              indy(1),indy(2),indy(3),indy(4),
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              indx1, x2(ix2), nx2a(1), nx2b(1),
     *              nx2a(2), nx2b(2), nx2a(3), nx2b(3), x2a(1),
     *              x2b(1),x2a(2),x2b(2), x2a(3),x2b(3),
     *              indx2(1),indx2(2),indx2(3), cind, 
     *              flopcount, scr(iscr1), scr(iscr2))

         elseif (n1 .eq. 6 .and. n2 .eq. 4 .and. ncind .eq. 3) then
c           return 
            call cwork4643(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3),nya(4),nyb(4), 
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),indy(1),indy(2),indy(3),indy(4),
  
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              nx1a(5),nx1b(5),nx1a(6),nx1b(6),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4),x1a(5),x1b(5),x1a(6),x1b(6), 
     *              indx1, 
  
     *              x2(ix2), nx2a(1), nx2b(1),
     *              nx2a(2), nx2b(2), 
     *              nx2a(3), nx2b(3), 
     *              nx2a(4), nx2b(4), 
     *              x2a(1), x2b(1),x2a(2),x2b(2), 
     *              x2a(3), x2b(3),x2a(4),x2b(4), 
  
     *              indx2, cind, 
     *              flopcount, scr(iscr1), scr(iscr2), scr(iscr3))

         elseif (n1 .eq. 6 .and. n2 .eq. 2 .and. ncind .eq. 2) then
c           return 
            call cwork4622(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3),nya(4),nyb(4), 
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),indy(1),indy(2),indy(3),indy(4),
  
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              nx1a(5),nx1b(5),nx1a(6),nx1b(6),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4),x1a(5),x1b(5),x1a(6),x1b(6), 
     *              indx1, 
  
     *              x2(ix2), nx2a(1), nx2b(1),
     *              nx2a(2), nx2b(2), 
     *              x2a(1), x2b(1),x2a(2),x2b(2), 
  
     *              indx2, cind, 
     *              flopcount, scr(iscr1), scr(iscr2), scr(iscr3))


c-------------------------------------------------------------------------
c   Contracting into a 4-index output array.
c-------------------------------------------------------------------------

         elseif (n1 .eq. 4 .and. n2 .eq. 2 .and. ncind .eq. 1) then
            if (trace .and. 
     *          and(tracelevel, contraction_trace) .ne. 0) 
     *         print *,'cwork4421 at instruction ',current_line

c VFL 
c Make sure that the first index of  x2 is contracted.  


            call cwork4421(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3),nya(4),nyb(4), 
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),indy(1),indy(2),indy(3),indy(4),
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4), indx1, x2(ix2), nx2a(1), nx2b(1),
     *              nx2a(2), nx2b(2), x2a(1),
     *              x2b(1),x2a(2),x2b(2), indx2(1),indx2(2), cind, 
     *              flopcount, scr(iscr1), scr(iscr2))

c VFL 

         elseif (n1 .eq. 2 .and. n2 .eq. 4 .and. ncind .eq. 1) then
            if (trace .and. 
     *          and(tracelevel, contraction_trace) .ne. 0) 
     *         print *,'cwork4421 at instruction ',current_line

            call cwork4421(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3),nya(4),nyb(4), 
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),indy(1),indy(2),indy(3),indy(4),
     *              x2(ix2), nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *              nx2a(3),nx2b(3),nx2a(4),nx2b(4),
     *              x2a(1),x2b(1),x2a(2),x2b(2),x2a(3),x2b(3),
     *              x2a(4),x2b(4), indx2, 
     *              x1(ix1), nx1a(1), nx1b(1),
     *              nx1a(2), nx1b(2), x1a(1),
     *              x1b(1),x1a(2),x1b(2), indx1(1),indx1(2), cind, 
     *              flopcount, scr(iscr1), scr(iscr2))

         else if (n1 .eq. 4 .and. n2 .eq. 4 .and. ncind .eq. 2) then
             if (trace .and. 
     *           and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork4442 at instruction ',current_line

c--------------------------------------------------------------------------
c   The new cwork4442 always assumes the 1st index of "y" shows up in the
c   first operand array.  Therefore we must call it with the proper 
c   sequence of operands.
c--------------------------------------------------------------------------

            if (iseq .eq. 1) then
               call cwork4442(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),indy(1),indy(2),indy(3),indy(4),
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4), indx1, x2(ix2), 
     *              nx2a(1),nx2b(1),nx2a(2),nx2b(2),nx2a(3),nx2b(3),
     *              nx2a(4),nx2b(4),    
     *              x2a(1),
     *              x2b(1),x2a(2),x2b(2), x2a(3),x2b(3),x2a(4),x2b(4),
     *              indx2, cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))
            else
               call cwork4442(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),indy(1),indy(2),indy(3),indy(4),
     *              x2(ix2), nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *              nx2a(3),nx2b(3),nx2a(4),nx2b(4),
     *              x2a(1),x2b(1),x2a(2),x2b(2),x2a(3),x2b(3),
     *              x2a(4),x2b(4), indx2, x1(ix1),
     *              nx1a(1),nx1b(1),nx1a(2),nx1b(2),nx1a(3),nx1b(3),
     *              nx1a(4),nx1b(4),
     *              x1a(1),
     *              x1b(1),x1a(2),x1b(2), x1a(3),x1b(3),x1a(4),x1b(4),
     *              indx1, cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))
            endif
         else
            print *,'Error: Invalid contraction'
            print *,'ny, n1, n2, ncind = ',ny, n1, n2, ncind
            call abort_job()
         endif
      else if (ny .eq. 6) then
         if (n1 .eq. 4 .and. n2 .eq. 4 .and. ncind .eq. 1) then
               call cwork6441(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4), nya(5),nyb(5),
     *              nya(6),nyb(6),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),ya(5),yb(5),ya(6),yb(6),
     *              indy(1),indy(2),indy(3),indy(4),indy(5),indy(6),
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4), indx1, x2(ix2), 
     *              nx2a(1),nx2b(1),nx2a(2),nx2b(2),nx2a(3),nx2b(3),
     *              nx2a(4),nx2b(4),    
     *              x2a(1),
     *              x2b(1),x2a(2),x2b(2), x2a(3),x2b(3),x2a(4),x2b(4),
     *              indx2, cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))

         elseif (n1 .eq. 6 .and. n2 .eq. 4 .and. ncind .eq. 2) then
               call cwork6642(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4), nya(5),nyb(5),
     *              nya(6),nyb(6),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),ya(5),yb(5),ya(6),yb(6),
     *              indy(1),indy(2),indy(3),indy(4),indy(5),indy(6),
c 
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              nx1a(5),nx1b(5),nx1a(6),nx1b(6),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4), x1a(5),x1b(5),
     *              x1a(6),x1b(6), indx1, 
c 
     *              x2(ix2), 
     *              nx2a(1),nx2b(1),nx2a(2),nx2b(2),nx2a(3),nx2b(3),
     *              nx2a(4),nx2b(4),    
     *              x2a(1),x2b(1),x2a(2),x2b(2),x2a(3),x2b(3),
     *              x2a(4),x2b(4),
     *              indx2, cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))

         elseif (n1 .eq. 4 .and. n2 .eq. 6 .and. ncind .eq. 2) then
c         return 
               call cwork6642(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4), nya(5),nyb(5),
     *              nya(6),nyb(6),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),ya(5),yb(5),ya(6),yb(6),
     *              indy(1),indy(2),indy(3),indy(4),indy(5),indy(6),
  
     *              x2(ix2), nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *              nx2a(3),nx2b(3),nx2a(4),nx2b(4),
     *              nx2a(5),nx2b(5),nx2a(6),nx2b(6),
     *              x2a(1),x2b(1),x2a(2),x2b(2),x2a(3),x2b(3),
     *              x2a(4),x2b(4), x2a(5),x2b(5),
     *              x2a(6),x2b(6), indx2, 
  
     *              x1(ix1), 
     *              nx1a(1),nx1b(1),nx1a(2),nx1b(2),nx1a(3),nx1b(3),
     *              nx1a(4),nx1b(4),    
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4),
     *              indx1, 
  
     *              cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))

         elseif (n1 .eq. 6 .and. n2 .eq. 2 .and. ncind .eq. 1) then
               call cwork6621(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4), nya(5),nyb(5),
     *              nya(6),nyb(6),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),ya(5),yb(5),ya(6),yb(6),
     *              indy(1),indy(2),indy(3),indy(4),indy(5),indy(6),
c 
     *              x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *              nx1a(5),nx1b(5),nx1a(6),nx1b(6),
     *              x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *              x1a(4),x1b(4), x1a(5),x1b(5),
     *              x1a(6),x1b(6), indx1, 
c 
     *              x2(ix2), 
     *              nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *              x2a(1),x2b(1),x2a(2),x2b(2),
     *              indx2, cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))

         elseif (n1 .eq. 2 .and. n2 .eq. 6 .and. ncind .eq. 1) then
               call cwork6621(y(iy),nya(1),nyb(1), nya(2),nyb(2),
     *              nya(3),nyb(3), nya(4),nyb(4), nya(5),nyb(5),
     *              nya(6),nyb(6),
     *              ya(1),yb(1), ya(2),yb(2),ya(3),yb(3),
     *              ya(4),yb(4),ya(5),yb(5),ya(6),yb(6),
     *              indy(1),indy(2),indy(3),indy(4),indy(5),indy(6),
c 
     *              x2(ix2), nx2a(1),nx2b(1),nx2a(2),nx2b(2),
     *              nx2a(3),nx2b(3),nx2a(4),nx2b(4),
     *              nx2a(5),nx2b(5),nx2a(6),nx2b(6),
     *              x2a(1),x2b(1),x2a(2),x2b(2),x2a(3),x2b(3),
     *              x2a(4),x2b(4), x2a(5),x2b(5),
     *              x2a(6),x2b(6), indx2, 
c 
     *              x1(ix1), 
     *              nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *              x1a(1),x1b(1),x1a(2),x1b(2),
     *              indx1, cind, flopcount, scr(iscr1), scr(iscr2),
     *              scr(iscr3))


         else
            print *,'Error: Invalid contraction'
            print *,'ny, n1, n2, ncind = ',ny, n1, n2, ncind
            call abort_job()
         endif
      else if (ytype .eq. scalar_value) then

c--------------------------------------------------------------------------
c   Contraction of all indices into a scalar result.
c   1. x1 and x2 must have the same number of indices.
c   2. All indices of x1, x2 must be contracted out (i.e., are present in 
c      array "cind".
c--------------------------------------------------------------------------

         if (n1 .ne. n2 .or. n1 .ne. ncind) then
            print *,'Error: Invalid contraction'
            print *,'   Indices must match for contraction ',
     *                  'into scalar'
            print *,'   Contracted arrays are ',x1handle,x2handle
            print *,'   Arrays have ',n1,n2,' indices'
            print *,'   Num. of contracted indices ',ncind
            call abort_job()
         endif 

         iy = array_table(c_scalar_index, yhandle)
         if (n1 .eq. 2) then
            if (trace .and. 
     *          and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork0222 at instruction ',current_line

            call cwork0222(scalar_table(iy), x1(ix1),
     *                  nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *                  x1a(1),x1b(1),
     *                  x1a(2),x1b(2),
     *                  indx1,
     *                  x2(ix2),nx2a(1), nx2b(1), nx2a(2), nx2b(2),
     *                  x2a(1), x2b(1), x2a(2), x2b(2),
     *                  indx2(1), indx2(2), flopcount)
         else if (n1 .eq. 4) then
            if (trace .and.
     *          and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork0444 at instruction ',current_line
            call cwork0444(scalar_table(iy), 
     *           x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *           nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *           x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *           x1a(4),x1b(4), indx1, x2(ix2),
     *           nx2a(1),nx2b(1),nx2a(2),nx2b(2),nx2a(3),nx2b(3),
     *           nx2a(4),nx2b(4),
     *           x2a(1),
     *           x2b(1),x2a(2),x2b(2), x2a(3),x2b(3),x2a(4),x2b(4),
     *           indx2, flopcount)
         else if (n1 .eq. 3) then
            if (trace .and.
     *          and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork0333 at instruction ',current_line
            call cwork0333(scalar_table(iy),
     *           x1(ix1),nx1a(1),nx1b(1),nx1a(2),
     *           nx1b(2),nx1a(3),nx1b(3),
     *           x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *           indx1, x2(ix2),
     *           nx2a(1),nx2b(1),nx2a(2),nx2b(2),nx2a(3),nx2b(3),
     *           x2a(1),x2b(1),x2a(2),x2b(2), x2a(3),x2b(3),
     *           indx2, flopcount)
         else if (n1 .eq. 1) then
            if (trace .and.
     *          and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork0111 at instruction ',current_line
            call cwork0111(scalar_table(iy),
     *           x1(ix1),nx1a(1),nx1b(1),
     *           x1a(1),x1b(1),indx1, x2(ix2),
     *           nx2a(1),nx2b(1),x2a(1),x2b(1), 
     *           indx2, flopcount)
         else if (n1 .eq. 6) then
            if (trace .and.
     *          and(tracelevel, contraction_trace) .ne. 0)
     *          print *,'cwork0666 at instruction ',current_line
            call cwork0666(scalar_table(iy), 
     *           x1(ix1), nx1a(1),nx1b(1),nx1a(2),nx1b(2),
     *           nx1a(3),nx1b(3),nx1a(4),nx1b(4),
     *           nx1a(5),nx1b(5),nx1a(6),nx1b(6),
     *           x1a(1),x1b(1),x1a(2),x1b(2),x1a(3),x1b(3),
     *           x1a(4),x1b(4),x1a(5),x1b(5),x1a(6),x1b(6),indx1, 
     *           x2(ix2),
     *           nx2a(1),nx2b(1),nx2a(2),nx2b(2),nx2a(3),nx2b(3),
     *           nx2a(4),nx2b(4),nx2a(5),nx2b(5),nx2a(6),nx2b(6),
     *           x2a(1),
     *           x2b(1),x2a(2),x2b(2), x2a(3),x2b(3),x2a(4),x2b(4),
     *           x2a(5),x2b(5),x2a(6),x2b(6), 
     *           indx2, flopcount)
         else
            print *,'Error: Invalid scalar contraction'
            print *,'   n1, n2 = ',n1,n2
            call abort_job()
         endif

      else
         print *,'Error: Output of contraction must have 4 or less ',
     *           'indices, or be a scalar.'
         call abort_job()
      endif

      return
      end
