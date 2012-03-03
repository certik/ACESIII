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
      integer function find_current_block_map(x, array_table_entry,
     *                                  index_table, nindex_table,
     *                                  block_map_table, 
     *                                  nblock_map_table)
c-----------------------------------------------------------------------------
c   Looks up the current segments of each index of an array, then uses
c   those segments to find the block with those segments in the 
c   block_map_table.  
c-----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'

      integer x, nindex_table, nblock_map_table
      integer array_table_entry(larray_table_entry)
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)

      integer nindex, i, j, nb, iblock_map, nmatch, maxblk
      integer index(mx_array_index), seg(mx_array_index)
      integer lookup, block_map_lookup

      iblock_map = array_table_entry(c_block_map)
      if (iblock_map .eq. 0) then   ! array is not mapped

c----------------------------------------------------------------------------
c   Get a dummy block number. 
c----------------------------------------------------------------------------

         blkmgr_next_block = blkmgr_next_block + 1
         find_current_block_map = blkmgr_next_block
         return
      endif

      nindex = array_table_entry(c_nindex)
      do i = 1, nindex
         index(i) = array_table_entry(c_index_array1+i-1)
         seg(i)   = index_table(c_current_seg,index(i))
      enddo

      lookup = block_map_lookup(seg, nindex, x,
     *                         array_table_entry,
     *                         index_table, nindex_table)

      find_current_block_map = lookup - iblock_map + 1
     
      return
      end

