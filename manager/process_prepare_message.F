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
      subroutine process_prepare_message(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a prepare message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer node
      integer j, iblock
      integer*8 indblk, get_index_from_base
      integer ptr, msgbuffer, state, ierr
      integer status(MPI_STATUS_SIZE)
      integer diskloc, ifile, size, istat
      integer find_free_diskloc
      integer f_backupram
      logical flag, done
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then

c---------------------------------------------------------------------------
c   Post a recv for the prepare data into the dedicated message buffer.
c---------------------------------------------------------------------------

         indblk = get_index_from_base(base_mem_addr, x, 2)
         msgbuffer = server_msg(c_msg_msgbuffer,node)

         indblk = indblk + (msgbuffer-1)*server_mem_blocksize

         call mpi_irecv(x(indblk), server_msg(c_msg_size,node),
     *          MPI_DOUBLE_PRECISION, server_msg(c_msg_source,node),
     *          server_msg(c_msg_tag,node), mpi_comm_world,
     *          server_msg(c_msg_request,node), ierr)
         server_msg(c_msg_state,node) = recv_block_state
      endif

      if (state .eq. recv_block_state) then
         if (server_msg(c_msg_request,node) .eq. 
     *         MPI_REQUEST_NULL) then
            flag = .true.
         else
            call mpi_test(server_msg(c_msg_request,node), flag, 
     *                 status, ierr)
	 endif

         if (flag) then
            server_msg(c_msg_state,node) = wait_for_block_state
            server_msg(c_msg_cause,node) = null_cause
         endif
      endif

      if (server_msg(c_msg_state,node) .eq. wait_for_block_state) then
         call claim_memory_block(node, server_table, nserver_table,
     *                           .false.)
         if (server_msg(c_msg_state,node) .eq. null_state) then

c--------------------------------------------------------------------------
c   Swap the pointers of the dedicated message buffer and the data block.
c--------------------------------------------------------------------------

            msgbuffer = server_msg(c_msg_msgbuffer,node)
            iblock    = server_msg(c_msg_memptr,node)
 
c----------------------------------------------------------------------------
c   Change the server_table_ptr entries.
c----------------------------------------------------------------------------

            ptr = server_table_ptr(iblock)
            if (ptr .gt. 0) then

c----------------------------------------------------------------------------
c   If this block is dirty, we are simply going to overwrite it.  No disk
c   backup is done.  Therefore, we must mark it as "clean".
c----------------------------------------------------------------------------

               if (and(server_table(c_server_flags,ptr),
     *             server_dirty_flag) .ne. 0)  then
                  call mark_block_clean(iblock, server_table, 
     *                               nserver_table)
               endif

               server_table(c_server_flags,ptr)  = 0  ! clear the flag
               server_table(c_server_memloc,ptr) = 0
               server_table(c_server_busy_node,ptr) = 0
            endif

c----------------------------------------------------------------------------
c   Swap the pointers within the message node.
c----------------------------------------------------------------------------

            server_msg(c_msg_msgbuffer,node) = iblock
            server_table_ptr(iblock)         = -1
            server_msg(c_msg_memptr,node)    = msgbuffer

c----------------------------------------------------------------------------
c   Find the server_table_ptr for the new memptr (i. e. the data we just 
c   received.  Mark the new block as dirty.
c----------------------------------------------------------------------------

            ptr = server_msg(c_msg_stptr,node)
            iblock = server_msg(c_msg_memptr,node)
            server_table_ptr(iblock)          = ptr
            server_table(c_server_memloc,ptr) = iblock
            call mark_block_dirty(iblock, server_table, nserver_table)

c---------------------------------------------------------------------------
c   Mark the message as complete and return.
c---------------------------------------------------------------------------

            server_msg(c_msg_state,node) = null_state
            server_msg(c_msg_cause,node) = null_cause
         endif
      endif 
      return
      end
