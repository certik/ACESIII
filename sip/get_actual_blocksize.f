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
      subroutine get_actual_blocksize(array, array_block, blkndx,
     *              array_table, narray_table,
     *              index_table, nindex_table,
     *              segment_table, nsegment_table, nwblock)
c--------------------------------------------------------------------------
c   Returns the actual size of data in a given block of an array.
c   No check is made to determine whether or not the block exists.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer array, array_block, narray_table, nindex_table,
     *        nsegment_table, nwblock, blkndx
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer nind, iblk, i, val1, val2
      integer ind(mx_array_index), seg(mx_array_index)

      nind = array_table(c_nindex, array)
      iblk = blkndx
      if (iblk .lt. 0) then
         print *,'Error: get_actual_blocksize could not find array ',
     *         array,' block ',array_block
         call dump_block_ids()
         call abort_job()
      endif
 
      call get_block_indices(iblk, ind)
      call get_block_segments(iblk, seg)

      nwblock = 1
      do i = 1, nind
         call get_index_segment(ind(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
         nwblock = nwblock * (val2-val1+1)
      enddo

      return
      end
