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
      integer function block_map_lookup(seg, nindex, array, 
     *                         array_table_entry,
     *                         index_table, nindex_table)
c-------------------------------------------------------------------------
c   Returns a lookup index into the block_map table for the block of a
c   given array that contains data for the segments (seg(i), i = 1, nindex).
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      integer nindex, array
      integer seg(nindex)
      integer array_table_entry(larray_table_entry)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
     
      integer i, n, index, eseg
      integer offset
      integer nseg(mx_array_index), bseg(mx_array_index)

c-------------------------------------------------------------------------
c   Look up the number of segments in each index.
c-------------------------------------------------------------------------

      do i = 1, nindex
         index = array_table_entry(c_index_original+i-1)
         bseg(i) = index_table(c_bseg, index)
         eseg    = index_table(c_eseg, index) 
         nseg(i) = eseg - bseg(i) + 1
      enddo

      offset = array_table_entry(c_block_map) + seg(1)-bseg(1)
      n = 1
      do i = 2, nindex
         n = n * nseg(i-1) 
         offset = offset + n * (seg(i) - bseg(i))
      enddo 

      block_map_lookup = offset
      return
      end
