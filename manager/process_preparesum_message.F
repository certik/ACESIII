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
      subroutine process_preparesum_message(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a prepare message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer node
      integer i, j, iblock
      integer*8 indblk, ind1, ind2, get_index_from_base
      integer ptr, iorequest, msgbuffer, state, ierr
      logical flag, done
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer status(mpi_status_size)
      integer*8 n64, ixx

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
         ind1= indblk + (msgbuffer-1)*server_mem_blocksize
c         print *,'Server ',me,' RECV: node ',node,' msgbuf ',
c     *      msgbuffer,' source ',server_msg(c_msg_source,node),
c     *     ' array ',
c     *     server_msg(c_msg_array,node),' segs ',
c     *    (server_msg(c_msg_bsegs+j-1,node),
c     *     server_msg(c_msg_esegs+j-1,node),j=1,4),
c     *     ' size ',server_msg(c_msg_size,node),' seqno ',
c     *     server_msg(c_msg_seqno,node),' tag ',
c     *     server_msg(c_msg_tag,node)

         call mpi_irecv(x(ind1), server_msg(c_msg_size,node),
     *          MPI_DOUBLE_PRECISION, 
     *          server_msg(c_msg_source,node),
     *          server_msg(c_msg_tag,node), mpi_comm_world,
     *          server_msg(c_msg_request,node), ierr)
         server_msg(c_msg_state,node) = recv_block_state
      endif

      if (server_msg(c_msg_state,node) .eq. recv_block_state) then

c----------------------------------------------------------------------------
c   Test for completion of the data transfer.
c----------------------------------------------------------------------------

         if (server_msg(c_msg_request,node) .eq. 
     *       MPI_REQUEST_NULL) then
            flag = .true.
         else
            call mpi_test(server_msg(c_msg_request,node), flag,
     *                 status, ierr)
         endif

         if (flag) then
c            print *,'Server ',me,' PREPARESUM RECV COMPLETE',
c     *         ': node ',node,
c     *         ' msgbuffer ',server_msg(c_msg_msgbuffer,node),
c     *         ' iblock ',server_msg(c_msg_memptr,node)
            server_msg(c_msg_state,node) = wait_for_block_state
         endif
      endif

      if (server_msg(c_msg_state,node) .eq. wait_for_block_state) then
         call claim_memory_block(node, server_table, nserver_table,
     *                           .true.)

         if (server_msg(c_msg_state,node) .eq. null_state) then

c---------------------------------------------------------------------------
c   All necessary data is now in memory.  Compute the sum if necessary.
c---------------------------------------------------------------------------

            iblock = server_msg(c_msg_memptr,node)
            msgbuffer = server_msg(c_msg_msgbuffer, node)

            if (server_msg(c_msg_flag,node) .eq. 
     *                       preparesum_copy_flag) then
               
c---------------------------------------------------------------------------
c   No need to sum, simply change the pointers to make the message buffer
c   available as the data buffer.
c---------------------------------------------------------------------------

               ptr = server_msg(c_msg_stptr,node)
               server_table_ptr(iblock) = -1  ! iblock becomes msg buffer
               server_msg(c_msg_msgbuffer,node) = iblock
               server_table_ptr(msgbuffer) = ptr
               server_msg(c_msg_memptr,node) = msgbuffer
               server_table(c_server_memloc,ptr) = msgbuffer

               iblock = server_msg(c_msg_memptr,node)
               indblk = get_index_from_base(base_mem_addr, x, 2)
               ind1 = indblk + (iblock-1) * server_mem_blocksize

c                print *,'Server ',me,' Swap data in node ',node,
c     *                ' new iblock ',
c     *                server_msg(c_msg_memptr,node),' contents ',
c     *                x(ind1), ' new msgbuffer ',
c     &                server_msg(c_msg_msgbuffer,node),
c     *                ' new pointer target for ptr ',ptr,' is ',
c     *                server_table(c_server_memloc,ptr),' diskloc ',
c     *                server_table(c_server_diskloc,ptr)
            else

c--------------------------------------------------------------------------
c   Sum the 2 buffers.
c--------------------------------------------------------------------------

               indblk = get_index_from_base(base_mem_addr, x, 2)
               ind1 = indblk + (iblock-1) * server_mem_blocksize
               ind2 = indblk + (msgbuffer-1)* server_mem_blocksize

               do i = 1, server_msg(c_msg_size,node)
                  x(ind1+i-1) = x(ind1+i-1) + x(ind2+i-1)
               enddo
            endif

c--------------------------------------------------------------------------
c   The block is now dirty, since we have put new data into it.
c--------------------------------------------------------------------------

            call mark_block_dirty(server_msg(c_msg_memptr,node), 
     *                            server_table, nserver_table)
         endif
      endif
      return
      end
