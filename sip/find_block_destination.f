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
      integer function find_block_destination(x, array_table, 
     *               narray_table, index_table, nindex_table, 
     *               block_map_table, nblock_map_table)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'

      integer x, narray_table, nindex_table, nblock_map_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,nblock_map_table)

      integer i, nindex, iblock_map, nb, type
      integer ii, j, nmatch
      integer index(mx_array_index), seg(mx_array_index)

      find_block_destination = -1

c-------------------------------------------------------------------------
c   Find the current segments of the array's indices.
c-------------------------------------------------------------------------

      type = array_table(c_array_type, x)
      if (type .ne. distributed_array .and.
     *    type .ne. served_array) then
         print *,'Error: find_block_destination was called with ',
     *     'an array that is not distributed or served.'
         print *,'Array ',x,' type = ',type
         call abort_job()
      endif

      nindex = array_table(c_nindex, x)
      do i = 1, nindex
        index(i) = array_table(c_index_array1+i-1,x)
        seg(i)   = index_table(c_current_seg,index(i))
      enddo

      iblock_map = array_table(c_block_map, x)
      nb         = array_table(c_numblks, x)
  
c----------------------------------------------------------------------------
c   Search the block map table for a block that matches the segments.
c----------------------------------------------------------------------------

      do i = 1, nb
         ii = iblock_map + i - 1
         nmatch = 0
         do j = 1, nindex
            if (seg(j) .eq. block_map_table(c_block_map_seg+j-1,ii))
     *          nmatch = nmatch + 1
         enddo

         if (nmatch .eq. nindex) then

c--------------------------------------------------------------------------
c   We have a match.  Return the processor.
c--------------------------------------------------------------------------

            find_block_destination = block_map_table(c_processor,ii)
            return
         endif
      enddo

      if (find_block_destination .eq. -1) then
         print *,'Error: Cannot find block destination ',
     *     ' for array ',x,' in block_map_table'
         print *,'Current index ',(index(i),i=1,nindex)
         print *,'Current segs: ',(seg(i),i=1,nindex)
         print *,'Current line ',current_line
         print *,'Block_map_table: nblock_map_table = ',
     *     nblock_map_table, ' nb = ',nb
         do i = 1, nb
            print *,(block_map_table(j,iblock_map+i-1),
     *            j=1,lblock_map_entry)
         enddo

         print *,'Matching seg = ',(seg(i),i=1,nindex)
         print *,'nmatch = ',nmatch
         call abort_job()
      endif
      return
      end
