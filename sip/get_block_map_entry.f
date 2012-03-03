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
      subroutine get_block_map_entry(array_table_entry, block, 
     *                               block_map_table, nblock_map_table,
     *                               block_map_entry)
      implicit none
      include 'interpreter.h'

      integer array_table_entry(larray_table_entry)
      integer nblock_map_table
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer block
      integer block_map_entry(lblock_map_entry)

      integer i, iblock_map

c------------------------------------------------------------------------
c   Get the block map entry and move it into the return argument.
c------------------------------------------------------------------------

      iblock_map = array_table_entry(c_block_map) + block - 1
      if (iblock_map .gt. nblock_map_table) then
         print *,'GET_BLOCK_MAP_ENTRY: Block_map_table has ',
     *            nblock_map_table,' but attempted reference for ',
     *            'block ',iblock_map
         print *,'Array_table pointer = ',
     *            array_table_entry(c_block_map),
     *           ' block in arg. list = ',block
         call abort_job()
      endif

      do i = 1, lblock_map_entry
         block_map_entry(i) = block_map_table(i,iblock_map)
      enddo
      return
      end
