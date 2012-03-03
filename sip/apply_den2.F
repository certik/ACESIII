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
      subroutine apply_den2(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Each element of the array source(i,j) is divided by the 
c   corresponding element of the array target(i,j).  
c   Format of command is :
c
c   execute apply_den2 source target
c
c   There is no resriction on the type of array of source or target. 
c
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      integer*8 indblk, get_block_index
      integer block, blkndx, seg
      integer find_current_block

      integer i, j, ii, source, target, company, comm, ierr
      integer source_type, target_type
      integer nindex_source, nindex_target, stack
      integer*8 isource, itarget 
      integer n(mx_array_index),ind(mx_array_index) 
      integer sval1(mx_array_index),sval2(mx_array_index)  
      integer tval1(mx_array_index),tval2(mx_array_index)  
      integer msg(len_sip_server_message)
      integer status(MPI_STATUS_SIZE) 

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif


c--------------------------------------------------------------------------
c   Locate the data for both arrays.  
c--------------------------------------------------------------------------

      target = op(c_op1_array)
      source = op(c_result_array)

c---------------------------------------------------------------------------
c   Look up source's address.
c---------------------------------------------------------------------------

      block  = find_current_block(source, array_table(1,source),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,source)
      isource = get_block_index(source, block, stack, blkndx, x, .true.) 

c---------------------------------------------------------------------------
c   Look up target's address.
c---------------------------------------------------------------------------

      block  = find_current_block(target, array_table(1,target),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,target)
      itarget = get_block_index(target, block, stack, blkndx, x, .true.) 

c--------------------------------------------------------------------------
c   Check the dimensions of both arrays both arrays.  
c--------------------------------------------------------------------------
      
      source_type = array_table(c_array_type, source)
      target_type = array_table(c_array_type, target)

      nindex_source = array_table(c_nindex, source)
      nindex_target = array_table(c_nindex, target)

      if (nindex_source .ne. nindex_target) then
         print *,'Error: Both arrays in apply_den2 must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex_source,' indices.'
         print *,'Second array has ',nindex_target,' indices.'
         call abort_job()
      endif

      if (nindex_source .ne. 2) then 
         print *,'Arrays in apply_den2 must have 2 indices' 
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source 
         ind(i) = array_table(c_index_array1+i-1,source)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, sval1(i), sval2(i))

      enddo

c-------------------------------------------------------------------------
c   Get segment ranges of the array target.
c-------------------------------------------------------------------------

      do i = 1, nindex_target 
         ind(i) = array_table(c_index_array1+i-1,target)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, tval1(i), tval2(i))

      enddo

      call do_div2(x(isource),x(itarget),
     *                      sval1(1),sval2(1),sval1(2),sval2(2),
     *                      tval1(1),tval2(1),tval1(2),tval2(2))

      return
      end

      subroutine do_div2(x,y,
     *                   a1,a2,b1,b2,i1,i2,j1,j2)
      implicit none 
      include 'interpreter.h'

      integer a1,a2,b1,b2
      integer i1,i2,j1,j2
      double precision x(a1:a2,b1:b2)
      double precision y(i1:i2,j1:j2)

      integer i,j,k,l 
      double precision val, div

      do i = i1, i2 
      do j = j1, j2 
         div = y(i,j) 
c        if (dabs(div) .gt. 1.0d-5) then 
            val = x(i,j)/div   
            x(i,j) = val 
c        endif 
      enddo 
      enddo 

      return 
      end 

