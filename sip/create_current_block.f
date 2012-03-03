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
      subroutine create_current_block(array,array_table, 
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op, 
     *                 release_flag, direct_flag, blk, ierr)
c---------------------------------------------------------------------------
c   Creates the current block of the array from scratch.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer array, ierr, narray_table, nindex_table, nsegment_table
      integer nblock_map_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer op(loptable_entry)

      integer allocate_block
      integer i, j, blk, stack, blkndx, nwblock
      integer find_current_block
      integer find_current_block_map
      integer nind, dummy
      integer ind(mx_array_index)
      integer iseg(mx_array_index)
      integer lookup, block_map_lookup
      logical release_flag
      logical direct_flag

      integer nb, offset 

      ierr = 0

c----------------------------------------------------------------------------
c   direct_flag = .true. implies that it is already known that the block
c   is not present and must be allocated.  Therefore we can skip the 
c   search for the block in this case.
c---------------------------------------------------------------------------

      if (.not. direct_flag) then
         blk = find_current_block(array, array_table(1,array),
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           block_map_table, blkndx)
         if (blk .gt. 0) then
            ierr = blkndx
            return   !  Block already exists.
         endif
      endif

c----------------------------------------------------------------------------
c   Block is not present.  Allocate a new block and set it up for this
c   operation.
c
c   Determine the size of the data to be allocated.
c---------------------------------------------------------------------------

      nind = array_table(c_nindex, array)
      do i = 1, nind
         ind(i) = array_table(c_index_array1+i-1, array)
      enddo

      call determine_current_block_size(ind, nind,
     *        index_table, nindex_table,
     *        segment_table, nsegment_table, nwblock)

c--------------------------------------------------------------------------
c   Block is not present.  Allocate it and compute it.
c--------------------------------------------------------------------------

      blk = find_current_block_map(array, array_table(1,array),
     *                   index_table, nindex_table,
     *                   block_map_table, nblock_map_table)
      
      ierr = allocate_block(array, blk, nwblock, array_table,
     *                      narray_table, index_table, nindex_table,
     *                      block_map_table)
      if (ierr .le. 0) then
         print *,'Create_current_block: ',
     *              'Cannot allocate block for array, block = ',
     *              array, blk,' release_flag = ',release_flag 
         call array_block_summary(array_table, narray_table)
         print *,'Current segments:'
         do i = 1, nindex_table
            print *,'   index ',i,' current segment ',
     *          index_table(c_current_seg,i)
         enddo
         call dump_block_ids()
         call abort_job()
      else
         blkndx = ierr
      endif

      call blkmgr_insert_block_in_list(array_table(c_block_list,array), 
     *              dummy, blkndx, c_block_list_ptr, .false.)

c--------------------------------------------------------------------------
c   Store the block indices in the block header.
c--------------------------------------------------------------------------

      call set_block_indices(array, blk, blkndx, array_table(1,array))

c--------------------------------------------------------------------------
c   Set the current segment settings in the block header.
c--------------------------------------------------------------------------

      call set_block_segments(array, blk, blkndx, index_table,
     *                           nindex_table)

c---------------------------------------------------------------------------
c   Zero out the block.
c---------------------------------------------------------------------------

      stack = array_table(c_array_stack,array)
c      call clear_block(array, blk, stack, blkndx, nwblock)

      call set_opblock(array, blk, blkndx, op)

c---------------------------------------------------------------------------
c   If this is a distributed array, store the blkndx in the block_map_table
c   entry.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,array) .eq. distributed_array) then 
         do i = 1, nind
            iseg(i) = index_table(c_current_seg,ind(i))
         enddo

         lookup = block_map_lookup(iseg, nind, array,
     *                          array_table(1,array), index_table, 
     *                          nindex_table)

         block_map_table(c_bmap_blkndx, lookup) = blkndx
      endif

      return
      end 

