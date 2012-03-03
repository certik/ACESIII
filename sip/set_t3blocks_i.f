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
      subroutine set_t3blocks_i(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c  The arrays to be used in the (T 4o) calculation are defined and their 
c  handles put into a common block.  
c  Input: array index of array to be used in compt3a  
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'saved_data.h'

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer array, array_type, nscalar_table, sind 
      double precision scalar_table(nscalar_table), tind 
      integer*8 address_table(narray_table)

      integer Ni, i, j, k
      integer flag, handle, index, nindex, ierr
      integer blk, blkndx

      integer t3iblocks  
      common /T3i_blocks/Ni, t3iblocks(20) 

c-------------------------------------------------------------------------
c   Put the array index in the appropriate slot.
c-------------------------------------------------------------------------

      array      = op(c_op1_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. scalar_value) then
         print *,'Error: The Second argument in set_t3_blocksa  
     *            must be a scalar.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in set_t3_blocksa, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      tind = scalar_table(sind)

      handle = op(c_result_array)
      if (tind .eq. 1.0) t3iblocks(1) = handle 
      if (tind .eq. 2.0) t3iblocks(2) = handle 
      if (tind .eq. 3.0) t3iblocks(3) = handle 
      if (tind .eq. 4.0) t3iblocks(4) = handle 
      if (tind .eq. 5.0) t3iblocks(5) = handle 
      if (tind .eq. 6.0) t3iblocks(6) = handle 
      if (tind .eq. 7.0) t3iblocks(7) = handle 
      if (tind .eq. 8.0) t3iblocks(8) = handle 
      if (tind .eq. 9.0) t3iblocks(9) = handle 
      if (tind .eq. 10.0) t3iblocks(10) = handle 
      if (tind .eq. 11.0) t3iblocks(11) = handle
      if (tind .eq. 12.0) t3iblocks(12) = handle
      if (tind .eq. 13.0) t3iblocks(13) = handle
      if (tind .eq. 14.0) t3iblocks(14) = handle
      if (tind .eq. 15.0) t3iblocks(15) = handle
      if (tind .eq. 16.0) t3iblocks(16) = handle
      if (tind .eq. 17.0) t3iblocks(17) = handle
      if (tind .eq. 18.0) t3iblocks(18) = handle
      if (tind .eq. 19.0) t3iblocks(19) = handle
      if (tind .eq. 20.0) t3iblocks(20) = handle

      if (tind .gt. 20.0) then 
         write(6,*) ' Only 20 arrays are currently allowed in 
     *                set_t3blocks_i' 
         write(6,*) ' There are at least :', tind 
         call abort_job()
      endif

      return
      end

