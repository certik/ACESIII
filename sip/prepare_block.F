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
      subroutine prepare_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer)
c---------------------------------------------------------------------------
c   Client-side implementation of the prepare instruction.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'proto_events.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_barrier_data.h'
      include 'server_monitor.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer instruction_timer, comm_timer

      integer i, j, k
      integer handle, op1, op1_type, ierr
      integer comm_flag, whandle
      integer my_rank, wincomm, my_group, mpi_group_world
      integer nxt
      integer blkndx
      
      integer msg(len_sip_server_message,mx_msg)
      integer imsg, find_free_msg_buffer

      integer*8 indblk, get_block_index
      integer stack
      integer blk
      integer find_current_block
      integer find_block_destination
      integer request
      integer ind(mx_array_index), seg(mx_array_index)
      integer status(MPI_STATUS_SIZE)
      integer nind, blocksize
      integer get_block_request, block_map_lookup
      integer mgr
      integer f_form_msg_tag, tag
      integer mx_msg_used
      character*2 optype
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
      save msg

      handle = op(c_result_array)
      op1    = op(c_op1_array)
      op1_type = array_table(c_array_type,op1)
      if (array_table(c_array_type,handle) .ne. served_array) then
         print *,'Error: prepare operation must use served array'
         print *,'Array = ',handle,' type = ',
     *            array_table(c_array_type,handle)
         print *,'Operation: ',(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      mx_msg_used = mx_msg
#ifdef ALTIX
      dptr = dshptr
c      mx_msg_used = min(mx_msg, 100)
#endif

c--------------------------------------------------------------------------
c   Set the "prepare" flag, indicating that the server_barrier should destroy
c   all local copies of the array's blocks.
c--------------------------------------------------------------------------

      array_table(c_prepare_flag,handle) = 1

c--------------------------------------------------------------------------
c   Get the segments of the array.
c--------------------------------------------------------------------------

      nind = array_table(c_nindex,handle)
      do i = 1, nind
         ind(i) = array_table(c_index_array1+i-1,handle)
         seg(i) = index_table(c_current_seg,ind(i))
      enddo

c---------------------------------------------------------------------------
c   Find the data to be sent to the I/O mgr.
c---------------------------------------------------------------------------

      blk = find_current_block(op1, array_table(1,op1),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)
      if (blk .lt. 0) then
         print *,'Task ',me,' Cannot find current block of',
     *      ' operand array in PREPARE instruction.'
         print *,'Array = ',op1,' current segments ',
     *         (seg(i),i=1,nind)
         call array_block_summary(array_table, narray_table)
         call dump_block_ids()
         call abort_job()
      endif

      mgr = find_block_destination(handle, array_table,
     *               narray_table, index_table, nindex_table,
     *               block_map_table, nblock_map_table)

      stack  = array_table(c_array_stack,op1)
      indblk = get_block_index(op1, blk, stack, blkndx, x, .true.)

c--------------------------------------------------------------------------
c   Build a prepare message.
c--------------------------------------------------------------------------

      call block_on_server_delay(mgr)
      imsg = find_free_msg_buffer(server_requests, mx_msg_used, 
     *              instruction_timer, comm_timer)

      if (op(c_opcode) .eq. prepare_op) then
         msg(1,imsg) = sip_server_prepare
         optype = 'p'
      else if (op(c_opcode) .eq. prepare_increment_op) then
         msg(1,imsg) = sip_server_prepare_increment
         optype = 'p+'
      else
         print *,'Error: Invalid opcode in prepare_block: ',op(c_opcode)
         call abort_job()
      endif
      msg(2,imsg) = handle
      msg(3,imsg) = nind
      tag    = f_form_msg_tag()
      msg(4,imsg) = tag
      msg(5,imsg) = current_line
      msg(6,imsg) = op(c_server_stat_key)
      nxt =  7
      do i = 1, mx_array_index
         if (i .le. nind) then
            msg(nxt,imsg) = ind(i)
         else
            msg(nxt,imsg) = 0
         endif
         nxt = nxt + 1
      enddo

      blocksize = 1
      do i = 1, nind
         call get_index_segment(ind(i), seg(i), segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, msg(nxt,imsg), 
     *                          msg(nxt+1,imsg))
         blocksize = blocksize * (msg(nxt+1,imsg)-msg(nxt,imsg)+1)
         nxt = nxt + 2
      enddo

      do i = nind+1, mx_array_index
         msg(nxt,imsg)   = 0
         msg(nxt+1,imsg) = 0
         nxt = nxt + 2
      enddo

      msg(7,imsg) = block_map_lookup(seg, nind, handle,
     *                         array_table(1,handle),
     *                         index_table, nindex_table)

      call mpi_isend(msg(1,imsg), len_sip_server_message, 
     *                 mpi_integer, mgr, sip_server_message,
     *                 mpi_comm_world, server_requests(imsg), ierr)

c----------------------------------------------------------------------------
c   Now send the data to be "prepared".
c----------------------------------------------------------------------------

      request = get_block_request(op1, blk, blkndx, request)

c----------------------------------------------------------------------------
c   If the block is engaged in communication from a previous request, 
c   we must wait for it to finish before continuing.
c----------------------------------------------------------------------------

      if (request .ne. MPI_REQUEST_NULL) then
         call wait_on_block(op1, blk, blkndx, op1_type, request,
     *             instruction_timer, comm_timer)
      endif

      call server_monitor_write_log(blkndx, optype, mgr)
      call mpi_isend(x(indblk), blocksize, mpi_double_precision,
     *                 mgr, tag,
     *                 mpi_comm_world, request, ierr) 

      call set_block_request(op1, blk, blkndx, request)
      if (request .ne. MPI_REQUEST_NULL) 
     *    call blkmgr_insert_block_in_list(comm_list_head, 
     *            comm_list_tail, blkndx, c_comm_list_ptr,
     *            .true.)

      return
      end
