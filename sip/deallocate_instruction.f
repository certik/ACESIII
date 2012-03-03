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
      subroutine deallocate_instruction(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, op, 
     *                      comm_timer, instruction_timer)
c---------------------------------------------------------------------------
c   Runtime implementation code for the DEALLOCATE instruction.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blkmgr.h'
      include 'checkpoint_data.h'

      integer narray_table, nindex_table, nsegment_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer block_map_table(lblock_map_entry,*)
      integer op(loptable_entry)
      integer comm_timer, instruction_timer

      integer array, type, next, flag, flagtest
      integer stack, i, j, array_data, blk, blkndx, request
      integer get_block_request
      integer*8 ixx
      integer*8 get_index_from_base
      integer bhdr(1) 

      array = op(c_result_array)
      type  = array_table(c_array_type,array)
      if (type .ne. local_array) then
         print *,'Error: Deallocate instruction requires a local_array'
         print *,'Array, type = ',array,type
         call abort_job()
      endif

      flagtest = or(block_busy_flag,block_persistence_flag)
      ixx      = get_index_from_base(addr_blk_header, bhdr, 1)
      stack    = array_table(c_array_stack, array)

c-----------------------------------------------------------------------------
c   Scan the linked list of blocks in the array, releasing the blocks as
c   we go.
c-----------------------------------------------------------------------------

      next = array_table(c_block_list, array)
 2000 continue
      if (next .ne. 0) then
         i = next
         array_data =
     *            bhdr(ixx+(c_array_handle-1)*blkmgr_blocks+i-1)
c         call get_blk_header(array_data, i, c_array_handle)

         if (array_data .eq. array) then

c--------------------------------------------------------------------------
c   The block matches the array.  We can get rid of it now.
c--------------------------------------------------------------------------

            call get_block_id(stack, i, array, blk)
            blkndx = i
            array_table(c_current_blkndx, array) = i  ! Update array table entry

            call clear_block_created_flag(array, blk, blkndx)

c--------------------------------------------------------------------------
c   Must not free this block if it is engaged in communication.
c--------------------------------------------------------------------------

            request = get_block_request(array, blk, blkndx)
            if (request .ne. MPI_REQUEST_NULL)
     *         call wait_on_block(array, blk, blkndx, type, request,
     *             instruction_timer, comm_timer)

c-------------------------------------------------------------------------
c   Find the next block in the linked list.
c-------------------------------------------------------------------------

            call get_blk_header(next, i,c_block_list_ptr)
            if (next .eq. i) then
               print *,'Error: Recursion in block list'
               print *,'array = ',array,' head = ',
     *             array_table(c_block_list,array)
               call abort_job()
            endif

            call_marker = 30303
            call free_block(array, blk, blkndx, array_table,
     *                      narray_table, index_table, nindex_table,
     *                      block_map_table)
         endif

         go to 2000
      endif

c---------------------------------------------------------------------------
c   We have traversed the list.  We have deallocated all the blocks of the 
c   array.
c
c   Now, remove the last allocate for this array from the checkpoint data.
c---------------------------------------------------------------------------

      do i = nactive_allocate_table, 1, -1
         if (active_allocate_table(i) .eq. array) then 
            
c--------------------------------------------------------------------------
c   Shift remaining table data over this entry. 
c--------------------------------------------------------------------------

            do j = i+1, nactive_allocate_table
               active_allocate_table(j-1) = active_allocate_table(j)
               active_allocate_op(j-1)    = active_allocate_op(j)
            enddo

            nactive_allocate_table = nactive_allocate_table - 1
            return
         endif
      enddo

      return
      end
