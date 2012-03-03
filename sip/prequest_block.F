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
      subroutine prequest_block(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table,nblock_map_table,
     *                      op, instruction_timer, comm_timer, 
     *                      direct_flag)
c----------------------------------------------------------------------------
c   Handles processing of a "prequest" instruction (partial request of a 
c   block of a served array).  The form of this instruction is
c      
c   prequest T(a,i,b,j) = S(a,p,b,q), where
c      T(a,i,b,j) contains the requested data to be transferred from server.
c      S(a,p,b,q) is the server-resident block from which to get the data.
c      a,p,b,q are declared indices of the served array S.
c      i, j are declared as "index" and are used to select the slice of
c           data to be transferred.     
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
      integer nind, nind2
      integer ierr, handle, source_array, status(MPI_STATUS_SIZE)
      integer response_req, request, get_block_request
      integer blk, blkndx, mgr
      integer blocksize
      integer*8 indblk, get_block_index
      integer stack
      integer find_block_destination
      integer ind(mx_array_index), seg(mx_array_index)
      integer ind2(mx_array_index), seg2(mx_array_index)
      integer saveseg(mx_array_index)

      integer msg(len_sip_server_prequest, mx_msg)
      integer imsg, find_free_msg_buffer
      integer mx_msg_used
      integer nseg2, jseg, source_type

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer nxt
      integer f_form_msg_tag2, tag
      integer flag

      save msg

      mx_msg_used = mx_msg
#ifdef ALTIX
      dptr = dshptr
c      mx_msg_used = min(mx_msg, 100)
#endif
      handle       = op(c_result_array)
      source_array = op(c_op1_array)
      source_type = array_table(c_array_type, source_array)
      if (source_type .ne. served_array) then
         print *,'Error: Source array in prequest must be a ',
     *           'served array.'
         call abort_job() 
      endif

c------------------------------------------------------------------------
c   Make sure the requested block exists.  create_current_block will
c   create the block if it does not exist, and simply return if the 
c   block is already present.
c------------------------------------------------------------------------

      call create_current_block(handle,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., direct_flag, blk, ierr)
      blkndx = ierr
 
c--------------------------------------------------------------------------
c   Check the block_request_outstanding_flag.
c   If a previous request instruction on this block has not been satisfied,
c   there is no need to process the current request.
c--------------------------------------------------------------------------

      call get_block_request_outstanding_flag(handle,blk, blkndx,
     *                                       flag)
      if (flag .ne. 0) return

c---------------------------------------------------------------------------
c   Check the block_persistence_flag.  If it is already turned on, the
c   block is available from a previous access, and there is no need to
c   request it again.
c---------------------------------------------------------------------------

      call get_block_persistence_flag(handle, blk, blkndx, flag)
      if (flag .ne. 0) then
         call set_opblock(handle, blk, blkndx, op)
         call set_block_computed_flag(handle, blk, blkndx, 1)
         return
      endif
      
c------------------------------------------------------------------------
c   Pick up the current indices and segments of the arrays.
c------------------------------------------------------------------------

      nind = array_table(c_nindex, handle)
      do i = 1, nind
         ind(i) = array_table(c_index_array1+i-1,handle)
         seg(i) = index_table(c_current_seg,ind(i))
      enddo

      nind2 = array_table(c_nindex, source_array)
      do i = 1, nind2
         ind2(i) = array_table(c_index_array1+i-1,source_array)
      enddo

      if (nind .ne. nind2) then
         print *,'Error in prequest: Arrays must have the same ',
     *           'number of indices.'
         print *,'Array 1 has ',nind,' indices, array 2 has ',
     *            nind2,' indices.'
         call abort_job()
      endif

c-----------------------------------------------------------------------
c  Find the non-matching indices of the 2 arrays.
c-----------------------------------------------------------------------

      do i = 1, nind
         saveseg(i) = index_table(c_current_seg,ind2(i)) ! save for later restoration.
         if (ind(i) .eq. ind2(i)) then
            seg2(i) = seg(i)
         else
            
c------------------------------------------------------------------------
c   Indices do not match.  Find the segment of the source array which 
c   contains the value of the segment for "handle".
c------------------------------------------------------------------------

            nseg2 = index_table(c_nsegments,ind2(i))
            do jseg = 1, nsegment_table
               if (segment_table(c_index,jseg) .eq. ind2(i)) then
                  k = jseg
                  go to 100
               endif
            enddo

  100       continue
            do j = k, k+nseg2-1
               if (segment_table(c_range1,j) .le. seg(i) .and.
     *             segment_table(c_range2,j) .ge. seg(i)) then
                  seg2(i) = segment_table(c_segment,j)
                  go to 200
               endif
            enddo 
  200       continue 
         endif
      enddo
 
c------------------------------------------------------------------------
c   Change the "current_seg" field of each index of the source array.
c------------------------------------------------------------------------

      do i = 1, nind2
         index_table(c_current_seg,ind2(i)) = seg2(i)
      enddo
      
c-------------------------------------------------------------------------
c   Find the served array block containing the indices of the target array.
c-------------------------------------------------------------------------

      mgr = find_block_destination(source_array,  
     *               array_table, narray_table, 
     *               index_table, nindex_table,
     *               block_map_table, nblock_map_table)

c--------------------------------------------------------------------------
c   Build and send the request message to the mgr.
c--------------------------------------------------------------------------

      imsg = find_free_msg_buffer(server_requests, mx_msg_used, 
     *                 instruction_timer, comm_timer)
      msg(1,imsg) = sip_server_prequest
      msg(2,imsg) = source_array
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
         call get_index_segment(ind(i), seg(i), segment_table,
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

      do i = 1, mx_array_index
         if (i .le. nind) then
            msg(nxt,imsg) = ind2(i)
         else
            msg(nxt,imsg) = 0
         endif
         nxt = nxt + 1
      enddo

      do i = 1, nind
         call get_index_segment(ind2(i), seg2(i), segment_table,
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

      call mpi_isend(msg(1,imsg), len_sip_server_prequest, mpi_integer,
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
      if (request .ne. MPI_REQUEST_NULL) 
     *   call wait_on_block(handle, blk, blkndx, served_array, request,
     *             instruction_timer, comm_timer)
 
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
            call set_opblock(handle, blk, blkndx, op)
            call set_block_computed_flag(handle, blk, blkndx, 1)

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

c--------------------------------------------------------------------------
c   Restore the original values of the source_array's indices.
c--------------------------------------------------------------------------

      do i = 1, nind2
         index_table(c_current_seg,ind2(i)) = saveseg(i)
      enddo

      return
      end
