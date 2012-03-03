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
      subroutine sip_server_work()
      implicit none
      include 'server.h'
      include 'parallel_info.h'
      integer*8 iserver_table, get_index_from_base 
      integer i
      integer itab(1)
#ifdef ALTIX
      include 'sheap.h'
    
      pointer (iptr, itab)
#endif

c---------------------------------------------------------------------------
c   Find the location of the server table in memory.
c---------------------------------------------------------------------------

#ifdef ALTIX
      iptr = ishptr
#endif
      iserver_table = get_index_from_base(server_table_base_addr, 
     *                                    itab, 1)

c-----------------------------------------------------------------------------
c   Enter the server work routine.
c-----------------------------------------------------------------------------

      call sip_server(itab(iserver_table), nserver_table_entries)

c------------------------------------------------------------------------------
c   Remove the server's files.
c------------------------------------------------------------------------------

      do i = 1, nserver_blocksizes
         call f_deletefile(server_filename(i), server_unit(i))
      enddo
      return
      end

      subroutine sip_server(server_table, nserver_table)
c--------------------------------------------------------------------------
c   Driver routine for the sip data server.
c--------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_ckpt_data.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'
      include 'server_delay.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table) 

      integer nactive, i, j, n, next, imsg, nind
      integer iblock, msgbuffer
      integer ierr, ptr
      integer ifile
      logical flag
      integer status(MPI_STATUS_SIZE)
      integer status2(MPI_STATUS_SIZE)
      integer msg(100)
      integer*8 msg8(50)
      equivalence (msg, msg8)

      integer*8 indblk, get_index_from_base
      integer lendata
      integer array, size, diskloc, nbackup_threshold
      integer source, tag
      integer flagval
      integer comm
      integer pst_get_company_comm
      integer find_server_table_ptr
      integer count
      integer str_trimlen

      double precision x(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, x)
      pointer (iptr, ibuf)
#else
      common x
      equivalence (ibuf, x)
#endif
      integer signal
      double precision tmsg, tmsg_sq, trestore

      save signal

#ifdef ALTIX
      dptr = dshptr
      iptr = ishptr

      server_msgs_before_delay = -1
#else
      server_msgs_before_delay = -1
#endif
      delay_flag = 0
      server_msg_count = 0

      nserved_arrays = 0
      do i = 1, nserver_table

c---------------------------------------------------------------------------
c   Create a table of served arrays.
c---------------------------------------------------------------------------

         array = server_table(c_server_array,i)
         do j = 1, nserved_arrays
            if (served_array_table(j) .eq. array) then
               served_numblocks(j) = served_numblocks(j) + 1
               go to 10
            endif
         enddo

c----------------------------------------------------------------------------
c   Create entries for a new array.
c---------------------------------------------------------------------------

         nserved_arrays = nserved_arrays + 1
         served_array_table(nserved_arrays)  = array
         served_array_status(nserved_arrays) = 0
         served_array_entry(nserved_arrays)  = i
         served_numblocks(nserved_arrays)    = 1
   10    continue

c---------------------------------------------------------------------------
c   Fill in the blocksize entry for each block in the server table.
c---------------------------------------------------------------------------

         nind = server_table(c_server_nind,i)
         n = 1
         do j = 1, nind
            n = n * (server_table(c_server_esegs+j-1,i) - 
     *               server_table(c_server_bsegs+j-1,i) + 1)
         enddo
         server_table(c_server_size,i) = n
         server_table(c_server_iblk,i) = -99999999

c---------------------------------------------------------------------------
c   Determine the appropriate file on which to store this block.
c---------------------------------------------------------------------------

         do j = 1, nserver_blocksizes
            if (server_blocksizes(j) .ge. n) then
               server_table(c_server_file,i) = j
               go to 50
            endif
         enddo
   50 continue
      enddo

      nactive = 0
      nbackup_threshold = nserver_memblocks / 10
  100 continue

      if (nactive .eq. 0 .and.
     *    delay_flag .eq. 1) then
         delay_flag = 0
         server_msg_count = 0
      endif

c---------------------------------------------------------------------------
c   Do we have a message buffer available?
c---------------------------------------------------------------------------

      if (nactive .lt. nmessage_buffers) then
         
c----------------------------------------------------------------------------
c   Probe for a "ready" message from any source.
c----------------------------------------------------------------------------

         call mpi_iprobe(MPI_ANY_SOURCE, server_readytag, 
     *                   MPI_COMM_WORLD, 
     *                   flag, status, ierr)

  500    continue
         if (flag) then
            source = status(MPI_SOURCE)
            call mpi_get_count(status, MPI_INTEGER, count, ierr)

c--------------------------------------------------------------------------
c  Receive the new message.
c--------------------------------------------------------------------------

            nactive = nactive + 1
            do i = 1, count
               msg(i) = 0
            enddo
    
            call mpi_recv(msg, count, MPI_INTEGER,
     *                    source,server_readytag,
     *                    MPI_COMM_WORLD, status2,ierr)
            if (msg(1) .lt. 0) then
               print *,'Server ',me,' Recieved invalid message: ',
     *              (msg(j),j=1,count)
               call server_abort_job(server_table, nserver_table) 
            endif

            server_seqno = server_seqno + 1 

c--------------------------------------------------------------------------
c   Set delay flag if necessary.
c--------------------------------------------------------------------------

            server_msg_count = server_msg_count + 1
            if (.not. barrier_in_progress .and. 
     *           delay_flag .eq. 0) then
               if (server_msg_count .eq. server_msgs_before_delay)
     *             delay_flag = 1
            endif

c--------------------------------------------------------------------------
c   Decode the message into a message node.
c-------------------------------------------------------------------------

            call pop_work_node(imsg)
            server_msg(c_msg_type,imsg) = msg(1)
            server_msg(c_msg_array,imsg) = msg(2)
            server_msg(c_msg_source,imsg)   = source
            server_msg(c_msg_nind,imsg)     = msg(3)
            server_msg(c_msg_tag,imsg)      = msg(4)
            server_msg(c_msg_seqno,imsg)    = server_seqno
            server_msg(c_msg_current_line,imsg) = msg(5)
            server_msg(c_msg_stat_key,imsg)     = msg(6)
            server_msg(c_msg_iblk,imsg)         = msg(7) 
            if (barrier_in_progress) barrier_seqno = server_seqno

            nind = server_msg(c_msg_nind,imsg)
            if (nind .lt. 1 .or. nind .gt. mx_array_index) then
               nind = 0
               server_msg(c_msg_stptr,imsg) = 0
               server_msg(c_msg_size,imsg)  = 0
            else
               next = 7+mx_array_index
               n = 1
               do i = 1, nind
                  server_msg(c_msg_bsegs+i-1,imsg) = msg(next)
                  server_msg(c_msg_esegs+i-1,imsg) = msg(next+1)
                  n = n * (msg(next+1)-msg(next)+1)
                  next = next + 2
               enddo

               server_msg(c_msg_size,imsg)   = n
            endif

            server_msg(c_msg_state,imsg)  = begin_state
            server_msg(c_msg_cause, imsg) = null_cause
            server_msg(c_msg_memptr,imsg) = 0
            server_msg(c_msg_flag,imsg)   = 0

            if (msg(1) .eq. server_blocks_to_list_msg .or. 
     *          msg(1) .eq. server_list_to_blocks_msg) then
               blk_to_list_offset(1,imsg) = msg8(2)
               blk_to_list_offset(2,imsg) = msg8(3)
               blk_to_list_offset(3,imsg) = msg8(4)
            endif

c------------------------------------------------------------------------
c   Decode a second set of indices/segment ranges if the messsage is 
c   a partial request (prequest).
c-------------------------------------------------------------------------

            if (msg(1) .eq. server_prequest_msg) then
               next = 7 + 4 * mx_array_index
               do i = 1, nind
                  server_msg(c_msg_bsegs2+i-1,imsg) = msg(next)
                  server_msg(c_msg_esegs2+i-1,imsg) = msg(next+1)
                  next = next + 2
               enddo
            endif

c--------------------------------------------------------------------------
c   Perform the server array checking to maintain data integrity.
c--------------------------------------------------------------------------

            call check_for_server_array(imsg, server_table, 
     *                                  nserver_table)

c--------------------------------------------------------------------------
c   Find the data block in the server_table.
c--------------------------------------------------------------------------

            if (nind .gt. 0) then
               ptr = find_server_table_ptr(imsg,server_table,
     *                                     nserver_table, .true.)
               server_msg(c_msg_stptr,imsg) = ptr
               server_table(c_server_iblk,ptr) = msg(7)   ! iblk entry
            endif

 
c---------------------------------------------------------------------------
c   Add to performance stats.
c---------------------------------------------------------------------------

            if (do_stats) then
               stat_key = server_msg(c_msg_stat_key,imsg)
               if (stat_key .gt. 0) then
                  if (server_msg(c_msg_type,imsg) .eq. 
     *                               server_prepare_msgtype) then    
                     sstat_nprep(stat_key) = sstat_nprep(stat_key) + 1
                  else if (server_msg(c_msg_type,imsg) .eq.
     *                               server_prepare_increment) then
                     sstat_nprepsum(stat_key) = sstat_nprepsum(stat_key) 
     *                                        + 1
                  else if (server_msg(c_msg_type,imsg) .eq.
     *                               server_request_msgtype) then
                     sstat_nreq(stat_key) = sstat_nreq(stat_key) + 1
                  else if (server_msg(c_msg_type,imsg) .eq.
     *                               server_prequest_msg) then
                     sstat_npreq(stat_key) = sstat_npreq(stat_key) + 1
                  endif

                  sstat_msg_time1(imsg) = mpi_wtime()
               endif
            endif

             call insert_work_node(imsg)   ! add to tail of work queue
         endif   ! flag
      endif

c----------------------------------------------------------------------------
c   Grab the head of the work queue entry.
c---------------------------------------------------------------------------

      call remove_work_head(imsg)
      if (imsg .le. 0) then
         if (barrier_in_progress .and.
     *       .not. flag .and.
     *       barrier_msg_count .eq. nbarrier_msgs) then

c---------------------------------------------------------------------------
c   Barrier exit logic.
c---------------------------------------------------------------------------

            if (.not. flag) then

c---------------------------------------------------------------------------
c   Sleep, then do 1 more probe just to make sure there are no more 
c   messages waiting.
c---------------------------------------------------------------------------

#ifndef __crayx1
               call f_usleep(1000)
               call mpi_iprobe(MPI_ANY_SOURCE, server_readytag,
     *                   MPI_COMM_WORLD,
     *                   flag, status, ierr)
               if (flag) go to 500
#endif
               
               barrier_in_progress = .false.
               barrier_seqno       = 0
               barrier_msg_count   = 0

c---------------------------------------------------------------------------
c   Clear delay flag settings.
c---------------------------------------------------------------------------

               delay_flag = 0
               server_msg_count = 0

c---------------------------------------------------------------------------
c   Clear all the array statuses.
c---------------------------------------------------------------------------

               do i = 1, nserved_arrays
                  served_array_status(i) = 0
               enddo

c--------------------------------------------------------------------------
c   Send an acknowledgement message back to the sender of the barrier.
c   The message tag is a server_barrier_signal to avoid confusion with 
c   other messages being sent back from the server.
c--------------------------------------------------------------------------

               signal = server_barrier_signal
               call mpi_send(signal, 1, MPI_INTEGER,
     *                 0, server_barrier_signal, mpi_comm_world, 
     *                 status, ierr)
               call prt_time('Server leaving barrier mode')
            endif
         else

c--------------------------------------------------------------------------
c   Clear delay flag if necessary.
c--------------------------------------------------------------------------

            if (delay_flag .eq. 1) then
               if (.not. barrier_in_progress .and. 
     *             .not. flag) then

c---------------------------------------------------------------------------
c   We are not in a barrier mode, all work in the message buffers has 
c   been cleared, and the probe indicates no messages waiting.  Thus, we
c   may clear the delay flag for this server.
c---------------------------------------------------------------------------

                  delay_flag = 0
                  server_msg_count = 0
               endif
            endif 
         endif   ! barrier_in_progress, etc...
      else

c----------------------------------------------------------------------------
c   Perform the next step in processing this data.
c----------------------------------------------------------------------------

         call process_work_queue_entry(imsg, server_table, 
     *                                  nserver_table)
         if (server_msg(c_msg_state,imsg) .eq. null_state) then

c---------------------------------------------------------------------------
c   Message is completed.  Capture performance data.
c---------------------------------------------------------------------------

            if (do_stats .and. stat_key .gt. 0) then
               tmsg = mpi_wtime() - sstat_msg_time1(imsg)

               if (server_msg(c_msg_type,imsg) .eq.
     *                               server_prepare_msgtype) then
                  sstat_tprep(stat_key) = sstat_tprep(stat_key) + 
     *                                     tmsg
                  sstat_tprep2(stat_key) = sstat_tprep2(stat_key) +
     *                                     tmsg_sq
               else if (server_msg(c_msg_type,imsg) .eq.
     *                               server_prepare_increment) then
                  sstat_tprepsum(stat_key) = 
     *                     sstat_tprepsum(stat_key) + tmsg
                  sstat_tprepsum2(stat_key) =
     *                     sstat_tprepsum2(stat_key) + tmsg_sq
               else if (server_msg(c_msg_type,imsg) .eq.
     *                               server_request_msgtype) then
                  sstat_treq(stat_key) = sstat_treq(stat_key) + tmsg
                  sstat_treq2(stat_key) = 
     *                   sstat_treq2(stat_key) + tmsg_sq
               else if (server_msg(c_msg_type,imsg) .eq.
     *                               server_prequest_msg) then
                  sstat_tpreq(stat_key) = sstat_tpreq(stat_key) + tmsg
                  sstat_tpreq2(stat_key) = 
     *                   sstat_tpreq2(stat_key) + tmsg_sq
               endif
            endif

            call push_work_node(imsg)   ! thru with this one.
            nactive = nactive - 1

c---------------------------------------------------------------------------
c   Reset the busy flag on both the data buffer (if any) and the msg buffer.
c---------------------------------------------------------------------------
               
               iblock = server_msg(c_msg_memptr,imsg)
               if (iblock .gt. 0) then
                  ptr = server_table_ptr(iblock)
                  if (ptr .gt. 0)  then
                    flagval = and(server_table(c_server_flags,ptr),
     *                         server_busy_flag)
                    if (flagval .ne. 0) 
     *                  server_table(c_server_flags,ptr) = 
     *                   xor(server_busy_flag,
     *                   server_table(c_server_flags,ptr))
                     server_table(c_server_busy_node,ptr) = 0
                  endif
               endif
     
               msgbuffer = server_msg(c_msg_msgbuffer,imsg)
               if (msgbuffer .gt. 0) then
                  ptr = server_table_ptr(msgbuffer)
                  if (ptr .gt. 0) then
                    flagval = and(server_table(c_server_flags,ptr),
     *                         server_busy_flag)
                    if (flagval .ne. 0)
     *                 server_table(c_server_flags,ptr) =
     *                 xor(server_busy_flag,
     *                 server_table(c_server_flags,ptr))
                     server_table(c_server_busy_node,ptr) = 0
                  endif
               endif 
         else if (server_msg(c_msg_state,imsg) .eq. quit_state) then
            go to 1000                  ! job is done
         else
        
c------------------------------------------------------------------------
c   Check the msg_flag field of the node to see if a restore
c   is required.
c------------------------------------------------------------------------

            if (server_msg(c_msg_cause,imsg) .eq. restore_cause) then 

c------------------------------------------------------------------------
c   Perform a restore operation on the "memptr" block.
c------------------------------------------------------------------------

               iblock = server_msg(c_msg_memptr,imsg)
               ptr = server_table_ptr(iblock)
               diskloc = server_table(c_server_diskloc,ptr)
               ifile   = server_table(c_server_file,ptr)

               indblk = get_index_from_base(base_mem_addr, x, 2)
               indblk = (iblock-1)*server_mem_blocksize + indblk
               size = server_table(c_server_size,ptr)

               if (do_stats .and. stat_key .gt. 0)
     *                trestore = mpi_wtime()
c               print *,'RESTORE block ',iblock
               call f_restoreram(server_unit(ifile), diskloc, 
     *                           server_blocksizes(ifile),
     *                           x(indblk), size)
               if (do_stats .and. stat_key .gt. 0) then
                  trestore = mpi_wtime() - trestore
                  sstat_trestore(stat_key) = 
     *                    sstat_trestore(stat_key) + trestore
                  sstat_trestore2(stat_key) =
     *              sstat_trestore2(stat_key) + trestore * trestore
                  sstat_nrestore(stat_key) = 
     *                  sstat_nrestore(stat_key) + 1 
               endif
               server_msg(c_msg_cause,imsg) = null_cause
            endif

c------------------------------------------------------------------------
c   Back up a block if necessary.
c   A backup is necessary if (1) message processing requires a clean
c   block and can't get it, or (2) the number of clean blocks drops below 
c   the numbr of idle message buffers.
c------------------------------------------------------------------------

             if (server_msg(c_msg_cause,imsg) .eq.
     *            backup_agent_cause ) then
                  iblock = server_msg(c_msg_memptr,imsg)
                  call do_backup(iblock, server_table, 
     *                              nserver_table)
c                  if (iblock .eq. -1) then
c                     call find_clean_block(iblock,server_table,
c     *                     nserver_table)
c                     if (iblock .eq. -1) then
c                        print *,'SERVER DEADLOCK'
c                        call server_abort_job(server_table,
c     *                                       nserver_table)
c                     else
c                        call mark_block_clean(iblock, server_table,
c     *                                   nserver_table)
c                     endif
c                  endif
                  server_msg(c_msg_cause,imsg) = null_cause
            endif

c-------------------------------------------------------------------------
c   Put the node back at the end of the work queue.
c-------------------------------------------------------------------------

            if (imsg .lt. 1 .or. imsg .gt. mx_server_msg) then
               print *,'Task ',me,' NODE OUT OF RANGE: ',imsg
               call server_abort_job(server_table, nserver_table) 
            endif
            call insert_work_node(imsg) ! go to the end of the line.
            call check_work_stack_for_delay(flag, delay_flag,
     *                  server_msg_count)
         endif
      endif 

      go to 100

 1000 continue

      if (do_stats) then
         
c---------------------------------------------------------------------------
c   Get the IOCOMPANY communicator.
c---------------------------------------------------------------------------

         comm = pst_get_company_comm(me)
         call mpi_comm_size(comm, my_company_size, ierr)
         call mpi_comm_rank(comm, my_company_rank, ierr)

c----------------------------------------------------------------------------
c   Perform reduction operations on server stats.
c----------------------------------------------------------------------------

         if (my_company_size .eq. 1) go to 2000

         call mpi_reduce(sstat_tprep, x, mx_stat_keys, 
     *                   MPI_DOUBLE_PRECISION, MPI_SUM, 
     *                   0, comm, ierr)
         
         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tprep(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tprep2, x, mx_stat_keys,
     *                   MPI_DOUBLE_PRECISION, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tprep2(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tprepsum, x, mx_stat_keys, 
     *                   MPI_DOUBLE_PRECISION, MPI_SUM, 
     *                   0, comm, ierr)
         
         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tprepsum(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tprepsum2, x, mx_stat_keys,
     *                   MPI_DOUBLE_PRECISION, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tprepsum2(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_treq, x, mx_stat_keys, 
     *                   MPI_DOUBLE_PRECISION, MPI_SUM, 
     *                   0, comm, ierr)
         
         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_treq(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_treq2, x, mx_stat_keys,
     *                   MPI_DOUBLE_PRECISION, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_treq2(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tpreq, x, mx_stat_keys, 
     *                   MPI_DOUBLE_PRECISION, MPI_SUM, 
     *                   0, comm, ierr)
         
         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tpreq(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tpreq2, x, mx_stat_keys,
     *                   MPI_DOUBLE_PRECISION, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tpreq2(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_trestore, x, mx_stat_keys, 
     *                   MPI_DOUBLE_PRECISION, MPI_SUM, 
     *                   0, comm, ierr)
         
         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_trestore(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_trestore2, x, mx_stat_keys,
     *                   MPI_DOUBLE_PRECISION, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_trestore2(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tbackup, x, mx_stat_keys, 
     *                   MPI_DOUBLE_PRECISION, MPI_SUM, 
     *                   0, comm, ierr)
         
         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tbackup(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_tbackup2, x, mx_stat_keys,
     *                   MPI_DOUBLE_PRECISION, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_tbackup2(i) = x(i)
            enddo
         endif

         call mpi_reduce(sstat_nprep, ibuf, mx_stat_keys,
     *                   MPI_INTEGER, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_nprep(i) = ibuf(i)
            enddo
         endif

         call mpi_reduce(sstat_nprepsum, ibuf, mx_stat_keys,
     *                   MPI_INTEGER, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_nprepsum(i) = ibuf(i)
            enddo
         endif

         call mpi_reduce(sstat_nreq, ibuf, mx_stat_keys,
     *                   MPI_INTEGER, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_nreq(i) = ibuf(i)
            enddo
         endif

         call mpi_reduce(sstat_npreq, ibuf, mx_stat_keys,
     *                   MPI_INTEGER, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_npreq(i) = ibuf(i)
            enddo
         endif

         call mpi_reduce(sstat_nrestore, ibuf, mx_stat_keys,
     *                   MPI_INTEGER, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_nrestore(i) = ibuf(i)
            enddo
         endif

         call mpi_reduce(sstat_nbackup, ibuf, mx_stat_keys,
     *                   MPI_INTEGER, MPI_SUM,
     *                   0, comm, ierr)

         if (my_company_rank .eq. 0) then
            do i = 1, mx_stat_keys
               sstat_nbackup(i) = ibuf(i)
            enddo
         endif 

c-------------------------------------------------------------------------
c   Send accumulated data to the master.
c-------------------------------------------------------------------------

 2000    continue
         if (my_company_rank .eq. 0) then
            lendata = 12*mx_stat_keys
            call mpi_send(sstat_tprep, lendata, MPI_DOUBLE_PRECISION,  
     *                  0, server_stat_data_tag, mpi_comm_world,
     *                  status, ierr) 
            lendata = 6*mx_stat_keys
            call mpi_send(sstat_nprep, lendata, MPI_INTEGER,
     *                  0, server_stat_data_tag, mpi_comm_world,
     *                  status, ierr)
         endif
      endif

c---------------------------------------------------------------------------
c   Remove checkpoint files if they exist.
c---------------------------------------------------------------------------

      n = str_trimlen(ckpt_dat_filename)
      inquire (file=ckpt_dat_filename,exist=flag)
      if (flag) call f_deletefile(ckpt_dat_filename(1:n)//char(0),-1)
      n = str_trimlen(ckpt_ndx_filename)
      inquire (file=ckpt_ndx_filename,exist=flag)
      if (flag) call f_deletefile(ckpt_ndx_filename(1:n)//char(0),-1)
      return
      end

