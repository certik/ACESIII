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
      subroutine block_map_table_setup(array_table, narray_table,
     *                                 index_table, nindex_table, 
     *                                 segment_table, nsegment_table,
     *                                 nblock_map_table)
c---------------------------------------------------------------------------
c   Scans the array table to determine the number of entries necessary 
c   for the block_map_table.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer narray_table,nindex_table,nsegment_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer nblock_map_table
     
      integer i, nind, k, nblks
      integer bseg, eseg, nseg, ind

      nblock_map_table = 0
      do i = 11, narray_table   ! 1st 10 arrays are pre-defined.
         if (array_table(c_array_type,i) .eq. distributed_array .or.
     *       array_table(c_array_type,i) .eq. served_array ) then

c-------------------------------------------------------------------------
c   Define the block_map_entry table for the array.
c-------------------------------------------------------------------------

            nind = array_table(c_nindex,i)
            nblks = 1
            do k = 1, nind
               ind = array_table(c_index_array1+k-1,i)
               bseg = index_table(c_bseg, ind)
               eseg = index_table(c_eseg, ind)
               nseg = eseg - bseg + 1
               nblks = nblks * nseg
            enddo

            nblock_map_table = nblock_map_table + nblks
         endif 
      enddo

      return
      end
