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
      subroutine process_request_message(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a request message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer node
      integer i, iblock
      integer*8 indblk, get_index_from_base
      integer*8 ind1, ind2
      integer ptr, ptr2, msgbuffer, state, ierr
      integer memloc, diskloc, size, ifile
      integer status(MPI_STATUS_SIZE)
      logical flag, restore_data_flag
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
      double precision trestore

#ifdef ALTIX
      dptr = dshptr
#endif

      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then

c---------------------------------------------------------------------------
c   Check the server_table to determine if the data is available on disk or
c   in memory.
c---------------------------------------------------------------------------

         ptr = server_msg(c_msg_stptr,node)
         memloc = server_table(c_server_memloc,ptr)
         diskloc = server_table(c_server_diskloc,ptr)

         if (memloc .eq. 0 .and. 
     *       diskloc .eq. 0) then

c--------------------------------------------------------------------------
c   Calculate the integrals, storing them in the message buffer..
c--------------------------------------------------------------------------

            indblk = get_index_from_base(base_mem_addr, x, 2)
            msgbuffer = server_msg(c_msg_msgbuffer,node)
            indblk = indblk + (msgbuffer-1)*server_mem_blocksize
c            call compute_integrals_wrapper(
c     *                          server_msg(c_msg_bsegs,node),
c     *                          server_msg(c_msg_esegs,node),
c     *                          server_msg(c_msg_nind,node),
c     *                          x(indblk))
            print *,'Task ',me,' REQUEST ERROR: Attempt to request a',
     *          ' non-existent block'
            print *,'node ',node,' array ',server_msg(c_msg_array,node)
            print *,'Requested at line ',
     *          server_msg(c_msg_current_line,node)
            call server_abort_job(server_table, nserver_table) 
            server_msg(c_msg_flag,node) = integral_calculation_flag

c---------------------------------------------------------------------------
c   Send the data to the source processor.
c--------------------------------------------------------------------------

            call mpi_isend(x(indblk), server_msg(c_msg_size,node),
     *            mpi_double_precision, server_msg(c_msg_source,node),
     *            server_msg(c_msg_tag,node),mpi_comm_world,
     *            server_msg(c_msg_request,node), ierr)
            server_msg(c_msg_state,node) = wait_for_send_state
        else if (memloc .eq. 0 .and. diskloc .ne. 0) then
           
c-------------------------------------------------------------------------
c   Read data from disk into the message buffer.
c-------------------------------------------------------------------------

           msgbuffer = server_msg(c_msg_msgbuffer,node)
           indblk    = get_index_from_base(base_mem_addr, x, 2)
           size      = server_table(c_server_size,ptr)
           ind2      = indblk +(msgbuffer-1)*server_mem_blocksize
           ifile     = server_table(c_server_file, ptr)
           diskloc   = server_table(c_server_diskloc,ptr)

           if (do_stats) trestore = mpi_wtime()
           call f_restoreram(server_unit(ifile), diskloc,
     *                           server_blocksizes(ifile),
     *                           x(ind2), size)
ccccc           server_msg(c_msg_flag,node) = restore_flag

           if (do_stats .and. stat_key .gt. 0) then
              trestore = mpi_wtime() - trestore
              sstat_trestore(stat_key) = sstat_trestore(stat_key) + 
     *                         trestore
              sstat_trestore2(stat_key) = sstat_trestore2(stat_key) + 
     *                         trestore * trestore
              sstat_nrestore(stat_key) = sstat_nrestore(stat_key) + 1
           endif

c---------------------------------------------------------------------------
c   Send the data to its target.
c---------------------------------------------------------------------------

            call mpi_isend(x(ind2), size, mpi_double_precision,
     *               server_msg(c_msg_source,node),
     *               server_msg(c_msg_tag,node),mpi_comm_world,
     *               server_msg(c_msg_request,node), ierr)
	   server_msg(c_msg_state,node) = wait_for_send_state
        else if (memloc .ne. 0) then
             server_msg(c_msg_state,node) = wait_for_block_state
             server_msg(c_msg_flag,node)  = 0
        endif
      endif    ! begin_state

      if (server_msg(c_msg_state,node) .eq. 
     *                           wait_for_block_state) then
 
c-------------------------------------------------------------------------
c   Check the flags to see if the data is free.  It could be busied by a
c   prepare or preparesum operation that is partially complete.
c-------------------------------------------------------------------------

c          call claim_memory_block(node, server_table, nserver_table,
c     *                            .true.)
c          if (server_msg(c_msg_state,node) .ne. null_state) return
c          server_msg(c_msg_flag,node) = 0

c--------------------------------------------------------------------------
c   The data is free and in memory now.  Find the pointer to the data block.
c--------------------------------------------------------------------------

ccc          iblock = server_msg(c_msg_memptr,node)
ccc          ptr    = server_table_ptr(iblock)
           ptr = server_msg(c_msg_stptr,node)
           iblock = server_table(c_server_memloc,ptr)

c----------------------------------------------------------------------------
c   Free the block (it may have been busied by the claim_memory_block 
c   subroutine.
c---------------------------------------------------------------------------

ccc          server_table(c_server_flags,ptr) = xor(server_busy_flag,
ccc     *              server_table(c_server_flags,ptr))
ccc          server_table(c_server_busy_node,ptr) = 0

c--------------------------------------------------------------------------
c   Copy the data to the message buffer.
c--------------------------------------------------------------------------

          msgbuffer = server_msg(c_msg_msgbuffer,node)
          indblk = get_index_from_base(base_mem_addr, x, 2)
          ind1   = indblk + (iblock-1)*server_mem_blocksize
          ind2   = indblk + (msgbuffer-1)*server_mem_blocksize
          size   = server_table(c_server_size,ptr)
        
          do i = 1, size
             x(ind2+i-1) = x(ind1+i-1)
          enddo
             
c---------------------------------------------------------------------------
c   Post a mpi_isend using the data buffer as the data address.
c---------------------------------------------------------------------------

         call mpi_isend(x(ind2), size, mpi_double_precision,
     *               server_msg(c_msg_source,node),
     *               server_msg(c_msg_tag,node),mpi_comm_world,
     *               server_msg(c_msg_request,node), ierr)
         server_msg(c_msg_state,node) = wait_for_send_state

	 if (size .ne. server_msg(c_msg_size,node)) then
            print *,'Server ',me,' Error: size > msg_size'
            print *,'ptr ',ptr,' array ',
     *        server_table(c_server_array,ptr),' segs ',
     *        (server_table(c_server_bsegs+i-1,ptr),
     *         server_table(c_server_esegs+i-1,ptr),i=1,4)
            print *,'msg array ',server_msg(c_msg_array,node),
     *       ' segs ',(server_msg(c_msg_bsegs+i-1,node),
     *          server_msg(c_msg_esegs+i-1,node),i=1,4)
            call server_abort_job(server_table,nserver_table)
         endif
      endif   ! wait_for_block_state

      if (server_msg(c_msg_state,node) .eq.
     *                         wait_for_send_state) then

c-------------------------------------------------------------------------
c   Test for completion of the send.
c-------------------------------------------------------------------------

         if (server_msg(c_msg_request,node) .eq. 
     *                               MPI_REQUEST_NULL) then
            flag = .true.
         else
            call mpi_test(server_msg(c_msg_request,node), flag,
     *                 status, ierr)
         endif

         if (flag) then
CXXXX            if (server_msg(c_msg_flag,node) .eq. 0) then
               server_msg(c_msg_state,node) = null_state   ! done
CXXXX            else
CXXXX               server_msg(c_msg_state,node) = request_cleanup_state
CXXXX            endif
         else 
            return 
         endif
      endif

      if (server_msg(c_msg_state,node) .eq. request_cleanup_state) then
         
c--------------------------------------------------------------------------
c   Look up the node's data in the server_table.  If its entry points to 
c   a memory block, then another message has brought this data into 
c   memory, and no further processing is required.
c--------------------------------------------------------------------------

         ptr = server_msg(c_msg_stptr,node)
         iblock = server_table(c_server_memloc, ptr)
         if (iblock .le. 0) then
            server_msg(c_msg_state,node) = null_state
            return
         endif

c--------------------------------------------------------------------------
c   Find a new buffer, and swap the message buffer and data buffers, so 
c   that the message buffer enters the processing flow as a data buffer.
c--------------------------------------------------------------------------

         call make_free_block(node, server_table, nserver_table)
         if (server_msg(c_msg_state,node) .eq. null_state) then
            iblock    = server_msg(c_msg_memptr,node)
            msgbuffer = server_msg(c_msg_msgbuffer,node)
            ptr2      = server_table_ptr(iblock)
            if (ptr2 .gt. 0) then
               server_table(c_server_memloc,ptr2) = 0
            endif

            server_table_ptr(iblock)          = -1        ! new msg buffer
            server_table_ptr(msgbuffer)       = ptr
            server_table(c_server_memloc,ptr) = msgbuffer
            server_msg(c_msg_msgbuffer,node)  = iblock
            server_msg(c_msg_memptr,node)     = msgbuffer
         endif
      endif

      return
      end
