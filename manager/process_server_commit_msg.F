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
      subroutine process_server_commit_msg(node, server_table,
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
      integer array, nindex, tag, f_form_msg_tag 
      integer msglen
      integer msg(100)

      integer ckpt_ndx, ckpt_dat
      integer nblk_file, header_entry, next
      integer*8 loc, index_file_diskloc, data_file_diskloc 
      integer*8 ndx_file_len, array_size
      integer str_trimlen
      logical created_file, match, flag

      integer index_file_entry_size
      parameter (index_file_entry_size = 3+mx_array_index*2)
      integer*8 index_file_entry(index_file_entry_size)

      double precision xxx,xsum

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

c--------------------------------------------------------------------------
c   Open the checkpoint index file.
c--------------------------------------------------------------------------

         n = str_trimlen(ckpt_ndx_filename)
         call f_openfile(ckpt_ndx_filename(1:n) // char(0), ckpt_ndx)
         if (ckpt_ndx .eq. -1) then
            print *,'Error: Server ',me,' failed to open checkpoint',
     *              ' index file during commit operation.'
            call server_abort_job(server_table, nserver_table)
         endif

c---------------------------------------------------------------------------
c   Write out the new ndx_file_header.  Once this is written to disk, the
c   disk pointers point to the last copy of the array written to disk, and
c   any previously written copies become free for reuse.  
c---------------------------------------------------------------------------

         loc = 0
         call f_write_disk(ckpt_ndx, loc, ndx_file_header, 
     *                     index_file_header_size) 

c--------------------------------------------------------------------------
c   Close the checkpoint file.
c--------------------------------------------------------------------------

         call f_close_file(ckpt_ndx)

c--------------------------------------------------------------------------
c   Update the free space structures.
c   Any free space candidates are now moved to the free space list, and 
c   the candidate list is cleared.
c--------------------------------------------------------------------------

         do i = 1, nfree_ndx_can
            do j = 1, nfree_space_ndx
               if (free_space_ndx(1,j) .eq. -1 .and.
     *             free_space_ndx(2,j) .eq. -1) then

c----------------------------------------------------------------------------
c   Use existing free space entry.
c----------------------------------------------------------------------------

                  free_space_ndx(1,j) = free_space_ndx_candidate(1,i)
                  free_space_ndx(2,j) = free_space_ndx_candidate(2,i)
                  go to 100
               endif 
            enddo
 
c----------------------------------------------------------------------------
c   Add to free space list.
c----------------------------------------------------------------------------

            nfree_space_ndx = nfree_space_ndx + 1
            if (nfree_space_ndx .gt. mx_ckpt_arrays) then
                print *,'Error: Free space entries are exhausted.'
                call server_abort_job(server_table, nserver_table)
            endif

            free_space_ndx(1,nfree_space_ndx) = 
     *                  free_space_ndx_candidate(1,i)
            free_space_ndx(2,nfree_space_ndx) = 
     *                  free_space_ndx_candidate(2,i)
  100       continue
         enddo

         nfree_ndx_can = 0

         do i = 1, nfree_dat_can
            do j = 1, nfree_space_dat
               if (free_space_dat(1,j) .eq. -1 .and.
     *             free_space_dat(2,j) .eq. -1) then

c----------------------------------------------------------------------------
c   Use existing free space entry.
c----------------------------------------------------------------------------

                  free_space_dat(1,j) = free_space_dat_candidate(1,i)
                  free_space_dat(2,j) = free_space_dat_candidate(2,i)
                  go to 200
               endif 
            enddo
 
c----------------------------------------------------------------------------
c   Add to free space list.
c----------------------------------------------------------------------------

            nfree_space_dat = nfree_space_dat + 1
            if (nfree_space_dat .gt. mx_ckpt_arrays) then
                print *,'Error: Free space data entries are exhausted.'
                call server_abort_job(server_table, nserver_table)
            endif

            free_space_dat(1,nfree_space_dat) = 
     *                  free_space_dat_candidate(1,i)
            free_space_dat(2,nfree_space_dat) = 
     *                  free_space_dat_candidate(2,i)
  200       continue
         enddo

         nfree_dat_can = 0

c--------------------------------------------------------------------------
c   Send an acknowledgement to the master.  This tells him that the 
c   server has finished his checkpoint operation.
c   The message tag for the ack. msg is server_checkpoint_msg.
c--------------------------------------------------------------------------

            call mpi_isend(msg, 1,
     *            mpi_integer, 0, server_commit_msg,
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

