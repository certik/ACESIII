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
      subroutine process_work_queue_entry(node, server_table, 
     *                                    nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of the server messages.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_stat.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      integer node
      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)

      integer msg_type, ierr
      integer i, array, ix, ptr, ptr2

      msg_type = server_msg(c_msg_type,node)

      if (msg_type .eq. server_request_msgtype) then

c-------------------------------------------------------------------------
c   REQUEST MESSAGE
c-------------------------------------------------------------------------

         call process_request_message(node, server_table, nserver_table)
      else  if (msg_type .eq. server_prepare_msgtype) then

c-------------------------------------------------------------------------
c   PREPARE MESSAGE
c--------------------------------------------------------------------------

         call process_prepare_message(node, server_table, nserver_table)
      else if (msg_type .eq. server_prepare_increment) then

c--------------------------------------------------------------------------
c   PREPARESUM MESSAGE
c--------------------------------------------------------------------------

         call process_preparesum_message(node, server_table, 
     *                       nserver_table)
      else if (msg_type .eq. server_prequest_msg) then

c-------------------------------------------------------------------------
c   PREQUEST MESSAGE
c-------------------------------------------------------------------------

         call process_prequest_message(node, server_table, 
     *                                 nserver_table)
      else if (msg_type .eq. server_quit_msgtype) then

c--------------------------------------------------------------------------
c   QUIT MESSAGE
c--------------------------------------------------------------------------

         server_msg(c_msg_state,node) = quit_state
      else if (msg_type .eq. server_barrier_signal) then

c--------------------------------------------------------------------------
c   BARRIER MESSAGE
c--------------------------------------------------------------------------

         barrier_in_progress = .true.
         barrier_seqno       = server_msg(c_msg_seqno,node)
         barrier_msg_count   = barrier_msg_count + 1

c--------------------------------------------------------------------------
c   Pull out the line number and translate it into a key for the server
c   data collection.
c--------------------------------------------------------------------------

         if (dbg) print *,'Server ',me,' Entering barrier logic, line ',
     *      server_msg(c_msg_tag,node), ' seqno ',
     *      barrier_seqno,' source ',server_msg(c_msg_source,node),
     *      ' barrier count = ',barrier_msg_count

         call prt_time('Server time')
         server_msg(c_msg_state,node) = null_state
      else if (msg_type .eq. server_copy_msg) then
c----------------------------------------------------------------------------
c   SERVER_COPY MESSAGE - server-side copy of an array to a new array.
c----------------------------------------------------------------------------
         call process_server_copy_message(node, server_table,
     *                 nserver_table)
      else if (msg_type .eq. server_delete_msg) then
c----------------------------------------------------------------------------
c   SERVER_DELETE MESSAGE - server-side delete of an array's blocks.  
c      The blocks are not actually removed, they simply become available
c      for use of other arrays. 
c----------------------------------------------------------------------------
         call process_server_delete_message(node, server_table,
     *                 nserver_table)
      else if (msg_type .eq. server_blocks_to_list_msg) then
         if (mpi_io_support) then
            call process_server_blocks_to_list_msg(node, server_table,
     *                 nserver_table)
         else
            call process_server_blocks_to_list_no_mpi_io(node, 
     *                 server_table, nserver_table)
         endif
      else if (msg_type .eq. server_list_to_blocks_msg) then
         if (mpi_io_support) then
            call process_server_list_to_blocks_msg(node, server_table,
     *                 nserver_table)
         else
            print *,'Error: Received read_list_to_blocks_msg, but ',
     *              'system does not support MPI_IO.'
            call server_abort_job(server_table, nserver_table)
         endif
      else if (msg_type .eq. server_checkpoint_msg) then
         call process_server_checkpoint_msg(node, server_table,
     *                 nserver_table)
      else if (msg_type .eq. server_restart_msg) then
         call process_server_restart_msg(node, server_table,
     *                 nserver_table)
      else if (msg_type .eq. server_commit_msg) then
         call process_server_commit_msg(node, server_table,
     *                 nserver_table)
      else
         print *,'Error: Invalid message type in server: ',
     *         msg_type
         call server_abort_job(server_table, nserver_table)
      endif
      
      return
      end
