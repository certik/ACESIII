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
      subroutine determine_current_block_size(ind, nind,
     *        index_table, nindex_table,
     *        segment_table, nsegment_table, size)
c----------------------------------------------------------------------------
c   Calculates the size of the block defined by the "current_seg" segments
c   in the index table.  The size of the block is returned in the "size"
c   argument.
c
c   Because this subroutine uses only information from the index_table and
c   segment_table, it may be used to calculate a block's size before it is
c   actually allocated.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      integer nind
      integer ind(nind)
      integer nindex_table, nsegment_table
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer size

      integer i, segment, val1, val2

      size = 1
      do i = 1, nind
         segment = index_table(c_current_seg,ind(i))
         call get_index_segment(ind(i), segment, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
         size = size * (val2-val1+1)
      enddo

      return
      end
