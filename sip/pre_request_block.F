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
      subroutine pre_request_block(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table,nblock_map_table,
     *                      op, instruction_timer, comm_timer, 
     *                      direct_flag)
c---------------------------------------------------------------------------
c   This routine iterates the indices used in the current array request.
c   Then if the indices are well-formed, it sets up a block for the data
c   that is being requested and issues a request to the server.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_barrier_data.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry,
     *                           nblock_map_table)
      integer instruction_timer, comm_timer
      logical direct_flag

      integer i, j, k
      integer nind
      integer ierr, handle, status(MPI_STATUS_SIZE)
      integer response_req, request, get_block_request
      integer blk, blkndx, mgr
      integer blocksize
      integer*8 indblk, get_block_index
      integer stack
      integer find_block_destination, block_map_lookup
      integer find_current_block, find_current_block_map
      integer allocate_block_light
      integer ind(mx_array_index), seg(mx_array_index)
      integer next_segs(mx_array_index)
      integer dummy

      integer msg(len_sip_server_message, mx_msg)
      integer imsg, find_free_msg_buffer
      integer mx_msg_used

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer nxt
      integer f_form_msg_tag2, tag
      integer flag
      integer save_op(2)
      integer nwblock
      logical return_flag
      logical is_prefetch_index

      save msg

      mx_msg_used = mx_msg
#ifdef ALTIX
      dptr = dshptr
c      mx_msg_used = min(mx_msg, 100)
#endif
      handle = op(c_result_array)
      if (array_table(c_array_type, handle) .ne. served_array) then
         print *,'Task ',me,' Error: Request instruction uses ',
     *           'an array that is not type "served".'
         print *,'Array ',handle,' type ',
     *            array_table(c_array_type, handle)
         call abort_job()
      endif

c------------------------------------------------------------------------
c   Pick up the current indices and segments of the array.
c------------------------------------------------------------------------

      nind = array_table(c_nindex, handle)
      return_flag = .false.
      do i = 1, nind
         ind(i) = array_table(c_index_array1+i-1,handle)
         seg(i) = index_table(c_current_seg,ind(i))
         if (is_prefetch_index(ind(i)) ) then
            next_segs(i) = index_table(c_next_seg,ind(i))
	    if (next_segs(i) .eq. 0) return_flag = .true.  ! no prefetch is possible.
            index_table(c_current_seg,ind(i)) = next_segs(i)
         else
            next_segs(i) = seg(i)   ! index is not prefetch, use current seg
         endif
      enddo

      if (return_flag) go to 1000

c------------------------------------------------------------------------
c   create_current_block will store the block info in the opblock and 
c   opblkndx fields of the instruction.  We don't want the prefetched 
c   block data stored, because this data has already been set by the 
c   REQUEST instruction.  So we save this data and restore it later.
c------------------------------------------------------------------------

       save_op(1) = op(c_opblock)
       save_op(2) = op(c_opblkndx)    

c------------------------------------------------------------------------
c   Make sure the requested block exists.  create_current_block will
c   create the block if it does not exist, and simply return if the 
c   block is already present.
c------------------------------------------------------------------------

      blk = find_current_block(handle, array_table(1,handle),
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           block_map_table, blkndx)

      if (blk .le. 0) then

c----------------------------------------------------------------------------
c   Block is not present.  Allocate a new block and set it up for this
c   operation.
c
c   Determine the size of the data to be allocated.
c---------------------------------------------------------------------------

         call determine_current_block_size(ind, nind,
     *        index_table, nindex_table,
     *        segment_table, nsegment_table, nwblock)

c--------------------------------------------------------------------------
c   Block is not present.  Allocate it and compute it.
c--------------------------------------------------------------------------

         blk = find_current_block_map(handle, array_table(1,handle),
     *                   index_table, nindex_table,
     *                   block_map_table, nblock_map_table)

         ierr = allocate_block_light(handle, blk, nwblock, array_table,
     *                      narray_table, index_table, nindex_table,
     *                      block_map_table)
         blkndx = ierr 
         if (ierr .le. 0) then

c--------------------------------------------------------------------------
c   No blocks available right now.  We will simply not do this pre_request.
c--------------------------------------------------------------------------

            go to 1000
         endif

c--------------------------------------------------------------------------
c   Finish setting up the block.
c--------------------------------------------------------------------------

         call blkmgr_insert_block_in_list(
     *              array_table(c_block_list,handle),
     *              dummy, blkndx, c_block_list_ptr, .false.)

c--------------------------------------------------------------------------
c   Store the block indices  and segments in the block header.
c--------------------------------------------------------------------------

         call set_block_indices(handle, blk, blkndx, 
     *                          array_table(1,handle))
         call set_block_segments(handle, blk, blkndx, index_table,
     *                           nindex_table)
         stack = array_table(c_array_stack,handle)
      endif

      mgr = find_block_destination(handle, array_table,
     *               narray_table, index_table, nindex_table,
     *               block_map_table, nblock_map_table)

c--------------------------------------------------------------------------
c   Check the block_request_outstanding_flag.
c   If a previous request instruction on this block has not been satisfied,
c   there is no need to process the current request.
c--------------------------------------------------------------------------

      call get_block_request_outstanding_flag(handle,blk, blkndx,
     *                                       flag)
      if (flag .ne. 0) then

c---------------------------------------------------------------------------
c   Check block_computed_flag: If it is not set, this was probably a block
c   requested from an earlier loop.  Therefore we must set the 
c   block_computed_flag here.
c---------------------------------------------------------------------------

         call get_block_computed_flag(handle, blk, blkndx, flag)
         if (flag .eq. 0) then
            call set_block_computed_flag(handle, blk, blkndx, 1)

c---------------------------------------------------------------------------
c   If the block_scrub_flag is set, this has been scrubbed from usage in
c   an earlier loop.  We must turn this flag off in order to prevent the
c   block from possibly being reused.
c---------------------------------------------------------------------------

            call set_block_scrub_flag(handle, blk, blkndx, 0)
         endif
      
         go to 1000
      endif

c---------------------------------------------------------------------------
c   Check the block_persistence_flag.  If it is already turned on, the
c   block is available from a previous access, and there is no need to
c   request it again.
c---------------------------------------------------------------------------

      call get_block_persistence_flag(handle, blk, blkndx, flag)
      if (flag .ne. 0) then
         call set_block_computed_flag(handle, blk, blkndx, 1)

c---------------------------------------------------------------------------
c   If the block_scrub_flag is set, this has been scrubbed from usage in
c   an earlier loop.  We must turn this flag off in order to prevent the 
c   block from possibly being reused.
c---------------------------------------------------------------------------

         call set_block_scrub_flag(handle, blk, blkndx, 0)
         go to 1000
      endif
      
c--------------------------------------------------------------------------
c   Build and send the request message to the mgr.
c--------------------------------------------------------------------------

      imsg = find_free_msg_buffer(server_requests, mx_msg_used, 
     *                 instruction_timer, comm_timer)
      msg(1,imsg) = sip_server_request
      msg(2,imsg) = handle
      msg(3,imsg) = nind
      tag    = f_form_msg_tag2()
      msg(4,imsg) = tag
      msg(5,imsg) = current_line
      msg(6,imsg) = op(c_server_stat_key)
      nxt = 7 
      do i = 1, mx_array_index
         if (i .le. nind) then
            msg(nxt,imsg) = ind(i)
         else
            msg(nxt,imsg) = 0
         endif
         nxt = nxt + 1
      enddo

      do i = 1, nind
         call get_index_segment(ind(i), next_segs(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, msg(nxt,imsg), 
     *                             msg(nxt+1,imsg))
         nxt = nxt + 2
      enddo

      do i = nind+1, mx_array_index
         msg(nxt,imsg)   = 0
         msg(nxt+1,imsg) = 0
         nxt = nxt + 2
      enddo

      msg(7, imsg) = block_map_lookup(next_segs, nind, handle,
     *                         array_table(1,handle),
     *                         index_table, nindex_table)

      call mpi_isend(msg(1,imsg), len_sip_server_message, mpi_integer,
     *               mgr, sip_server_message,
     *                 mpi_comm_world, server_requests(imsg), ierr)

c-------------------------------------------------------------------------
c   Post a recv for the block of data.
c   Check if the block is engaged in ongoing communication first.
c-------------------------------------------------------------------------

      call get_actual_blocksize(handle, blk, blkndx,
     *                  array_table, narray_table,
     *                  index_table, nindex_table, 
     *                  segment_table, nsegment_table, blocksize)
      stack  = array_table(c_array_stack,handle)
      indblk = get_block_index(handle, blk, stack, blkndx, x, .true.)

      request = get_block_request(handle, blk, blkndx)
      if (request .ne. MPI_REQUEST_NULL) then
         call wait_on_block(handle, blk, blkndx, served_array, request,
     *             instruction_timer, comm_timer)
      endif
 
      call mpi_irecv(x(indblk), blocksize, mpi_double_precision,
     *               mgr, tag, mpi_comm_world,
     *               request, ierr)
      call set_block_request(handle, blk, blkndx, request)
      if (request .ne. MPI_REQUEST_NULL) then
         call set_block_request_outstanding_flag(handle, blk, blkndx)
         call blkmgr_insert_block_in_list(comm_list_head, 
     *             comm_list_tail, blkndx, c_comm_list_ptr,
     *             .true.)
      endif

c----------------------------------------------------------------------------
c   Set the block_computed_flag.
c----------------------------------------------------------------------------

      call get_block_computed_flag(handle, blk, blkndx, flag)
      if (flag .eq. 0) then
         call get_block_created_flag(handle, blk, blkndx, flag)
         if (flag .eq. 0) then

c--------------------------------------------------------------------------
c   Set the block_persistence_flag.  This guarantees that the block will
c   remain available as long as possible before having to do an actual
c   new request from it's "home" server.
c
c   Setting the flag also updates its age, which is used in the allocation
c   scheme.
c--------------------------------------------------------------------------

            call set_block_persistence_flag(handle, blk, blkndx, 1)
            call blkmgr_insert_block_in_list(persistent_list_head,
     *          persistent_list_tail, blkndx, c_persistent_ptr,
     *          .true.)

c--------------------------------------------------------------------------
c   Set the block request flag for possible monitoring.
c--------------------------------------------------------------------------

            call set_block_request_flag(handle, blk, blkndx, 1)
         endif
      endif

 1000 continue

c---------------------------------------------------------------------------
c   Restore indices to original state.
c---------------------------------------------------------------------------

      do i = 1, nind
         index_table(c_current_seg,ind(i)) = seg(i)
      enddo

      return
      end
