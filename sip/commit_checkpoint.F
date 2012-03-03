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
      subroutine commit_checkpoint(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, ierr
      integer company_comm, pst_get_company_comm
      integer pst_get_company
      integer request(10000)
      integer statuses(MPI_STATUS_SIZE,10000)
      integer niocompany
      integer msg(5)
      logical search

      character*40 ckpt_msg

      if (.not. mpi_io_support) return  

c--------------------------------------------------------------------------
c   Convert the free space candidates to actual free space entries.
c--------------------------------------------------------------------------

      search = .true. 
      do i = 1, nfree_space_candidates

c--------------------------------------------------------------------------
c   Find a free space entry that has been marked for use.
c--------------------------------------------------------------------------

         if (search) then
            do j = 1, nfree_space
               if (free_space(j) .lt. 0) then
                  free_space(j) = free_space_candidates(i)
                  free_space_sizes(j) = candidate_sizes(i)
                  go to 100
               endif
            enddo
         endif
 
c---------------------------------------------------------------------------
c   No free space entry was found.  Add to the free space list.
c---------------------------------------------------------------------------

         nfree_space = nfree_space + 1
         free_space(nfree_space) = free_space_candidates(i)
         free_space_sizes(nfree_space) = candidate_sizes(i) 
         search = .false.
  100    continue
      enddo

      nfree_space_candidates = 0

      if (me .eq. 0) then

c--------------------------------------------------------------------------
c   Send commit signals to each server.
c--------------------------------------------------------------------------

         msg(1) = sip_server_commit_signal
         msg(4) = current_line   ! in the "tag" field of the msg.
         niocompany = 0
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Send the checkpoint message.
c---------------------------------------------------------------------------

               niocompany = niocompany + 1
               call mpi_isend(msg, 4,
     *            MPI_INTEGER, i-1,
     *            sip_server_message, mpi_comm_world,
     *            request(niocompany), ierr)
            endif
         enddo

c--------------------------------------------------------------------------
c   Write the master restart data to the master checkpoint file.
c--------------------------------------------------------------------------

         call write_master_checkpoint_data(scalar_table,
     *                   nscalar_table, index_table, nindex_table) 

c--------------------------------------------------------------------------
c   Wait on server requests to complete.
c--------------------------------------------------------------------------

         if (niocompany .gt. 0) then
            call mpi_waitall(niocompany, request, statuses, ierr)

c--------------------------------------------------------------------------
c   Post a recv for each server's acknowledgement signal.
c--------------------------------------------------------------------------

            niocompany = 0
            do i = 1, nprocs
               if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Receive the "checkpoint complete" message.
c   Each server will send a 1-word indicator message using the
c   sip_server_checkpoint_signal as the message tag.
c---------------------------------------------------------------------------

                  niocompany = niocompany + 1
                  call mpi_irecv(msg, 1,
     *               MPI_INTEGER, i-1,
     *               sip_server_commit_signal, mpi_comm_world,
     *               request(niocompany), ierr)
               endif
            enddo

c-----------------------------------------------------------------------
c   Wait for all the acknowledgement messages.
c------------------------------------------------------------------------

            call mpi_waitall(niocompany, request, statuses,
     *                    ierr)
         endif   ! niocompany .gt. 0
      endif   ! me .eq. 0

c--------------------------------------------------------------------------
c   Non-master members of the worker company must wait at this barrier
c   until the master arrives.
c--------------------------------------------------------------------------

      company_comm = pst_get_company_comm(me) 
      call mpi_barrier(company_comm, ierr)
      if (dbg) then
         print *,'Task ',me,' Passed checkpoint worker barrier'
         call c_flush_stdout()
      endif

      return
      end
