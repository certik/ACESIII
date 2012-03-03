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
      subroutine server_takedown(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op )
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
      integer request(nprocs-my_company_size)
      integer statuses(MPI_STATUS_SIZE,mx_msg)
      integer status2(MPI_STATUS_SIZE,nprocs-my_company_size)
      integer msg(5)

      if (io_company_id .eq. 0) return
      if (dbg) print *,'Task ',me,
     *   ' is waiting at Server takedown barrier: line number ',
     *   current_line

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
      endif

c--------------------------------------------------------------------------
c   Wait until the servers have received all the descriptor messages that
c   have been sent.
c--------------------------------------------------------------------------

      call mpi_waitall(mx_msg, server_requests, 
     *                    statuses, ierr)
      if (dbg) print *,'Task ',me,
     *     ' All server desc. messages have completed.' 

c--------------------------------------------------------------------------
c   The initial barrier guarantees that no more request/prepare instructions
c   will be initiated.
c---------------------------------------------------------------------------

      call mpi_barrier(company_comm, ierr)

c-------------------------------------------------------------------------
c   Send a quit message to each I/O server from rank 0 of this company.
c-------------------------------------------------------------------------

      if (my_company_rank .eq. 0) then
         niocompany = 0
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Send the barrier message.
c---------------------------------------------------------------------------

               niocompany = niocompany + 1
               call mpi_isend(sip_server_barrier_signal, 1, 
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
         if (dbg) print *,
     *     'All server barrier signals have been received.'

c-----------------------------------------------------------------------------
c   Now send each server a "quit" signal.
c-----------------------------------------------------------------------------

         niocompany = 0
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Send the quit message.
c---------------------------------------------------------------------------

               niocompany = niocompany + 1
               call mpi_isend(sip_server_quit, 1, 
     *               MPI_INTEGER, i-1, 
     *               sip_server_message, mpi_comm_world,
     *               request(niocompany), ierr)
            endif
         enddo

c---------------------------------------------------------------------------
c   Wait until all servers have received their quit signal.
c---------------------------------------------------------------------------

         call mpi_waitall(niocompany, request, 
     *                    statuses, ierr)
         if (dbg) print *,'All server quit signals have been received.'
      endif 

c--------------------------------------------------------------------------
c   Force all processors to await acknowledgement that the servers have
c   received the barrier signal.
c--------------------------------------------------------------------------

      call mpi_barrier(company_comm, ierr)
      if (dbg) print *,'Task ',me,' passed the server takedown barrier.'
      return
      end
