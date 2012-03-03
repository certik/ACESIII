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
      subroutine process_prequest_message(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a prequest message.
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
      double precision xmin, xmax

#ifdef ALTIX
      dptr = dshptr
#endif
      state = server_msg(c_msg_state,node)
c      print *,'PREQ MSG: line ',server_msg(c_msg_current_line,node),
c     *     ' node ',node,' state ',state,' cause ',
c     *     server_msg(c_msg_cause,node), 
c     *     ' clean ptr ',
c     *     clean_block_ptr,' nclean ',nclean_blocks
      if (state .eq. begin_state .or. 
     *    state .eq. wait_for_block_state) then

c---------------------------------------------------------------------------
c   Force the data into memory.  
c---------------------------------------------------------------------------
 
         call claim_memory_block(node, server_table, nserver_table,
     *                            .true.)

         if (server_msg(c_msg_state, node) .eq. null_state) then
            server_msg(c_msg_state, node) = prequest_intermediate_state
         else
c            print *,'CLAIM FAILED: cause ',server_msg(c_msg_cause,node)
            return   ! claim did not work, retry later.
         endif
      endif
          
      if (server_msg(c_msg_state,node) .eq. 
     *                 prequest_intermediate_state) then

c--------------------------------------------------------------------------
c   Lookup the addresses of the data block and the msgbuffer data block.
c--------------------------------------------------------------------------

         ptr       = server_msg(c_msg_stptr,node)
         memloc    = server_table(c_server_memloc,ptr)
         msgbuffer = server_msg(c_msg_msgbuffer,node)

         indblk    = get_index_from_base(base_mem_addr, x, 2)
         ind1      = indblk + (memloc-1)*server_mem_blocksize
         ind2      = indblk +(msgbuffer-1)*server_mem_blocksize

c---------------------------------------------------------------------------
c   Copy the slice of data from the data block (ind1) into the msgbuffer 
c   (ind2).
c---------------------------------------------------------------------------

         call prequest_copy_slice(node, x(ind1), x(ind2),
     *                   server_table, nserver_table)

c---------------------------------------------------------------------------
c   Post a mpi_isend using the msgbuffer as the data address.
c---------------------------------------------------------------------------

         size = server_msg(c_msg_size,node)
         call mpi_isend(x(ind2), size, mpi_double_precision,
     *               server_msg(c_msg_source,node),
     *               server_msg(c_msg_tag,node),mpi_comm_world,
     *               server_msg(c_msg_request,node), ierr)

c--------------------------------------------------------------------------
c   Release the data block by turning off the server_table flags.  This 
c   allows the block to be reused in other messages at this point.
c---------------------------------------------------------------------------

          server_table(c_server_flags,ptr) = xor(server_busy_flag,
     *              server_table(c_server_flags,ptr))
          server_table(c_server_busy_node,ptr) = 0
          call push_clean_block(memloc,server_table, nserver_table)

c---------------------------------------------------------------------------
c   Enter the "wait_for_send" state.
c---------------------------------------------------------------------------

         server_msg(c_msg_state,node) = wait_for_send_state
      endif

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
c            print *,'PREQ line ',server_msg(c_msg_current_line,node),
c     *            ' for node ',node,' completed!!!'
            server_msg(c_msg_state,node) = null_state   ! done
         endif
      endif

      return
      end
