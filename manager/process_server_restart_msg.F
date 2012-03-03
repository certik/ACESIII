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
      subroutine process_server_restart_msg(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_checkpoint_msg message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_ckpt_data.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer node
      integer i, j, n, istart, nsearch, iblock
      integer*8 indblk, get_index_from_base
      integer*8 ind, imsg
      integer msgbuffer, state, ierr
      integer memloc, diskloc, size, ifile, nxt
      integer request, status(MPI_STATUS_SIZE)
      integer array, iarray, nindex, tag, f_form_msg_tag 
      integer msglen
      integer msg(100)

      integer ckpt_ndx, ckpt_dat
      integer nblk_file, header_entry, next
      integer*8 loc, index_file_diskloc, data_file_diskloc 
      integer str_trimlen
      integer bsegs(mx_array_index), esegs(mx_array_index)
      logical created_file, match, flag

      integer index_file_entry_size
      parameter (index_file_entry_size = 3+mx_array_index*2)
      integer*8 index_file_entry(index_file_entry_size)
      integer ptr, find_server_table_ptr
      integer nind, iblk
      integer istat, find_free_diskloc, f_backupram


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
c   Determine a pointer to the server's message buffer.
c---------------------------------------------------------------------------

         ind = get_index_from_base(base_mem_addr, x, 2)
         msgbuffer = server_msg(c_msg_msgbuffer,node)

         ind = ind + (msgbuffer-1)*server_mem_blocksize

c---------------------------------------------------------------------------
c   Open the checkpoint files.
c----------------------------------------------------------------------------

         n = str_trimlen(ckpt_ndx_filename)
         call f_openfile(ckpt_ndx_filename(1:n) // char(0), ckpt_ndx)
         if (ckpt_ndx .eq. -1) then

c----------------------------------------------------------------------------
c   No server data was previously checkpointed.  Nothing to do.
c----------------------------------------------------------------------------

            go to 1000
         endif

         n = str_trimlen(ckpt_dat_filename)
         call f_openfile(ckpt_dat_filename(1:n) // char(0), ckpt_dat)
         if (ckpt_dat .eq. -1) then
            print *,'Error: Server ',me,' cannot open checkpoint ',
     *            'data file.'
            call abort_job()
         endif

c---------------------------------------------------------------------------
c   Read file header on the index file into memory.
c   For each entry in the file header:
c   Word 1 = array
c   Word 2 = number of blocks checkpointed by this server.
c   Word 3 = disk pointer to index file data.
c   Word 4 = number of indices in array.
c---------------------------------------------------------------------------

            loc = 0
            call f_read_disk(ckpt_ndx, loc, ndx_file_header,
     *                       index_file_header_size)

c----------------------------------------------------------------------------
c   Process each array in the index file header.
c----------------------------------------------------------------------------

            do iarray = 1, mx_ckpt_arrays
               array = ndx_file_header(1,iarray)
               if (array .ne. 0) then
                  nblk_file          = ndx_file_header(2,iarray)
                  index_file_diskloc = ndx_file_header(3,iarray) 
                  nind               = ndx_file_header(4,iarray)

                  do iblock = 1, nblk_file

c-------------------------------------------------------------------------
c   Read the index file entry for the block.
c-------------------------------------------------------------------------

                     call f_read_disk(ckpt_ndx, index_file_diskloc,
     *                        index_file_entry,
     *                        index_file_entry_size)
          
                     index_file_diskloc = index_file_diskloc + 
     *                                    index_file_entry_size 

c---------------------------------------------------------------------------
c   Retrieve the data from the index_file_entry array.
c---------------------------------------------------------------------------

                     iblk = index_file_entry(1)
                     size = index_file_entry(2)
                     data_file_diskloc = index_file_entry(3)
                     next = 4
                     do i = 1, mx_array_index
                       bsegs(i) =  index_file_entry(next) 
                       next = next + 1
                    enddo

                    do i = 1, mx_array_index
                       esegs(i) = index_file_entry(next) 
                       next = next + 1
                    enddo

c----------------------------------------------------------------------------
c   Create dummy info in the message node for use in the table search routine.
c----------------------------------------------------------------------------

                    server_msg(c_msg_type,node) =  server_restart_msg
                    server_msg(c_msg_array,node) = array
                    server_msg(c_msg_source,node) = 0
                    server_msg(c_msg_nind,node)   = nind
                    server_msg(c_msg_iblk,node)   = iblk
                    server_msg(c_msg_size,node)   = size
                    server_msg(c_msg_state,node)  = begin_state
                    server_msg(c_msg_cause,node)  = null_cause
                    server_msg(c_msg_memptr,node) = 0
                    server_msg(c_msg_flag,node)   = 0   
                    do i = 1, mx_array_index
                       if (i .gt. nind) then
                          server_msg(c_msg_bsegs+i-1,node) = 0
                          server_msg(c_msg_esegs+i-1,node) = 0
                       else
                          server_msg(c_msg_bsegs+i-1,node) = bsegs(i)
                          server_msg(c_msg_esegs+i-1,node) = esegs(i)
                       endif
                    enddo
  
c----------------------------------------------------------------------------
c   Find a match for this data in the server_table.
c----------------------------------------------------------------------------

                     ptr = find_server_table_ptr(node,server_table,
     *                                     nserver_table, .true.)
                     server_msg(c_msg_stptr,node) = ptr
                     server_table(c_server_iblk,ptr) = iblk

c---------------------------------------------------------------------------
c   Read the checkpoint data block from disk into the server's msg buffer.
c---------------------------------------------------------------------------

                     call f_read_disk(ckpt_dat, data_file_diskloc,
     *                        x(ind), size)
                     data_file_diskloc = data_file_diskloc + size

c---------------------------------------------------------------------------
c   Write the data to the server's SCR file.
c---------------------------------------------------------------------------

                     diskloc = server_table(c_server_diskloc,ptr)
                     ifile   = server_table(c_server_file,ptr)
                     if (diskloc .eq. 0) then
                        diskloc = find_free_diskloc(ifile, server_table,
     *                                   nserver_table) 
                        if (diskloc .lt. 0) then
                           next_server_diskloc(ifile) =
     *                          next_server_diskloc(ifile) + 1
                           diskloc = next_server_diskloc(ifile)
                        endif

                        server_table(c_server_diskloc,ptr) = diskloc
                     endif

                     istat =  f_backupram(server_unit(ifile), diskloc,
     *                                    server_blocksizes(ifile),
     *                                    x(ind), size)
                     if (istat .ne. 0)
     *                  call server_abort_job(server_table, 
     *                                        nserver_table)

                  enddo   ! iblock
               endif   ! array .ne. 0
            enddo   ! iarray
     
c--------------------------------------------------------------------------
c   Send an acknowledgement to the master.  This tells him that the 
c   server has finished his restart operation.
c   The message tag for the ack. msg is server_restart_msg.
c--------------------------------------------------------------------------

1000        continue
            call mpi_isend(msg, 1,
     *            mpi_integer, 0, server_restart_msg,
     *            mpi_comm_world,  server_msg(c_msg_request,node), 
     *            ierr)
            server_msg(c_msg_state,node) = wait_for_send_state
            return
      else if (state .eq. wait_for_send_state) then

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
               server_msg(c_msg_state,node) = null_state   ! done
         else
            return
         endif
      endif   ! state
      return
      end
