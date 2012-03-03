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
      subroutine server_barrier(array_table, narray_table,
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
      include 'server_barrier_data.h'
      include 'dbugcom.h'

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

      integer i
      integer ierr
      integer company_comm, pst_get_company_comm
      integer pst_get_company
      integer niocompany
      integer request(10000)
      integer msg(5)
      integer statuses(MPI_STATUS_SIZE,10000)
      integer status(MPI_STATUS_SIZE)

      if (dbg) then
         print *,'Task ',me,
     *    ' is waiting at Server barrier: line number ',
     *   current_line
         call prt_time('Worker time')
         call c_flush_stdout()
      endif

      company_comm = pst_get_company_comm(me)
      if (my_company_size .ne. 1) then 

c---------------------------------------------------------------------------
c   First, we do a sip_barrier.  This guarantees that no more request/prepare
c   instructions will be submitted from the workers.
c---------------------------------------------------------------------------

         call sip_barrier(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
      endif   ! my_company_size .ne. 1

c--------------------------------------------------------------------------
c   Wait until the servers have received all the descriptor messages that
c   have been sent.
c--------------------------------------------------------------------------

      call mpi_waitall(mx_msg, server_requests, 
     *                    statuses, ierr)
      if (dbg) print *,'Task ',me,' Completed waitall at line ',
     *     current_line 

c-------------------------------------------------------------------------
c   Send a barrier message to each I/O server from each worker of this company.
c-------------------------------------------------------------------------

      niocompany = 0
      msg(1) = sip_server_barrier_signal
      msg(4) = current_line   ! in the "tag" field of the msg.
      do i = 1, nprocs
         if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Send the barrier message.
c---------------------------------------------------------------------------

            niocompany = niocompany + 1
            call mpi_isend(msg, 4, 
     *               MPI_INTEGER, i-1, 
     *               sip_server_message, mpi_comm_world,
     *               request(niocompany), ierr)
         endif
      enddo

c---------------------------------------------------------------------------
c   Wait until all servers have received their barrier signal.
c---------------------------------------------------------------------------

      call mpi_waitall(niocompany, request, 
     *                    statuses, ierr)
      call prt_time('All server barrier signals were received.')

c---------------------------------------------------------------------------
c   Wait for a response from each server acknowledging receipt.
c---------------------------------------------------------------------------

      if (my_company_rank .eq. 0) then
         niocompany = 0
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Receive the barrier message.
c   Each server will send a 1-word indicator message using the 
c   sip_server_barrier_signal as the message tag.
c---------------------------------------------------------------------------

               niocompany = niocompany + 1
               call mpi_irecv(msg, 1,
     *               MPI_INTEGER, i-1,
     *               sip_server_barrier_signal, mpi_comm_world,
     *               request(niocompany), ierr)
            endif
         enddo

c-----------------------------------------------------------------------
c   Wait for all the acknowledgement messages.
c------------------------------------------------------------------------

         call mpi_waitall(niocompany, request, statuses, 
     *                    ierr)
      endif 

c---------------------------------------------------------------------------
c   Reset all persistent blocks  of all arrays which have had a "PREPARE"
c   since the last barrier.
c   This is done to insure data integrity.
c---------------------------------------------------------------------------

      do i = 1, narray_table
         if (array_table(c_prepare_flag,i) .ne. 0) then
            array_table(c_prepare_flag,i) = 0
            call free_persistent_blocks(i, array_table, narray_table,
     *               index_table, nindex_table, block_map_table)
         endif
      enddo

c--------------------------------------------------------------------------
c   Force all processors to await acknowledgement that the servers have
c   received the barrier signal.
c--------------------------------------------------------------------------

      call mpi_barrier(company_comm, ierr)
      call prt_time('Passed the server barrier.')
      return
      end
