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
      subroutine block_end_of_loop(array, block, blkndx,
     *             array_table, narray_table,
     *             index_table, nindex_table, block_map_table)
c--------------------------------------------------------------------------
c   Handles setting of flags at the end of a DO or PARDO loop.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blkmgr.h'

      integer array, block, blkndx
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer create_flag, flag
      integer get_block_request

c-------------------------------------------------------------------------
c   Check that this block is not the result of a CREATE instruction.  If it
c   is, it can only be released by DELETE'ing it.
c--------------------------------------------------------------------------

      call get_block_created_flag(array, block, blkndx, create_flag)
      if (create_flag .eq. 0) then
         call get_block_persistence_flag(array, block, blkndx,
     *                                   flag)
         call clear_block_computed_flag(array, block, blkndx)
         call set_block_busy_flag(array, block, blkndx, 0)

c---------------------------------------------------------------------------
c   This block is not a "persistent" block.
c---------------------------------------------------------------------------
            
         if (get_block_request(array, block, blkndx) .ne. 
     *                      MPI_REQUEST_NULL) then

c---------------------------------------------------------------------------
c   Mark the block as "scrubbable".  It will be freed later after its
c   communication is complete.
c---------------------------------------------------------------------------

            call set_block_scrub_flag(array, block, blkndx, 1)
c            call get_blk_header(flag, blkndx, c_block_flags)
c            print 10101,me,current_line,blkndx,flag
c10101 format('Task ',me,' line ',current_line,' blkndx, flag ',
c     *           2(i5,1x))
         else

c-------------------------------------------------------------------------
c   Free the block (if it is NOT a persistent block).
c-------------------------------------------------------------------------

            call_marker = 90909
            if (flag .eq. 0) 
     *         call free_block(array, block, blkndx, array_table,
     *                  narray_table, index_table, nindex_table,
     *                  block_map_table)
         endif
      endif

      return
      end
