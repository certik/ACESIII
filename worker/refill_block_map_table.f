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
      subroutine  refill_block_map_table(ncompany_workers,
     *                    block_map_table,
     *                    nblock_map_table,
     *                    array_table, narray_table,
     *                    index_table, nindex_table)
c---------------------------------------------------------------------------
c   Redistributes the distributed blocks of the block_map_table
c   as if the job is running with "ncompany_workers" processors.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      integer ncompany_workers
      integer narray_table, nindex_table, nblock_map_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)

      integer nblocks
      integer iarray
      integer i, j
      integer next
      integer next_worker
      
c--------------------------------------------------------------------------
c   Loop over each possible array.
c--------------------------------------------------------------------------

      next_worker = 0
      do iarray = 1, narray_table

c--------------------------------------------------------------------------
c   Search for a distributed array.
c--------------------------------------------------------------------------

         if (array_table(c_array_type,iarray) .eq. 
     *                                distributed_array) then

c---------------------------------------------------------------------------
c   Find the pointer to the first block of this array in the 
c   block_map_table.
c---------------------------------------------------------------------------

            next    = array_table(c_block_map,iarray)
            nblocks = array_table(c_numblks,iarray)

c--------------------------------------------------------------------------
c   Loop over each block of the array, filling in the block's new "home"
c   processor in the table.
c--------------------------------------------------------------------------

            do i = 1, nblocks
               block_map_table(c_processor,next) =
     *                          mod(next_worker,ncompany_workers)
               next        = next + 1
               next_worker = next_worker + 1
            enddo   ! nblocks
         endif      ! distributed array

c---------------------------------------------------------------------------
c   Modify entries for served array for new configuration.
c---------------------------------------------------------------------------

         if (array_table(c_array_type, iarray) .eq. 
     *                                     served_array) then
            nblocks = array_table(c_numblks,iarray)
            next    = array_table(c_block_map,iarray)

            do i = 1, nblocks
               block_map_table(c_processor,next) = ncompany_workers+1
               next = next + 1
            enddo
         endif
      enddo         ! iarray

      return
      end

