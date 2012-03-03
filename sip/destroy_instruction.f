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
      subroutine destroy_instruction(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer)
c--------------------------------------------------------------------------
c   Broadcast a server_delete message to all servers in the IOCOMPANY.
c   Format of command is :
c
c   destroy array
c
c   Should be executed after a server_barrier to insure that all messages 
c   referring to this array have been processed by all servers.  Otherwise
c   inconsistency of data may result.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer instruction_timer, comm_timer

      integer i, j, ii, array, company, comm, ierr
      integer array_type
      integer temp
      integer msg(len_sip_server_message)
      integer status(MPI_STATUS_SIZE) 
      integer pst_get_company
      integer pst_get_company_comm

      array = op(c_result_array)
      comm = pst_get_company_comm(me)

c--------------------------------------------------------------------------
c   The array must be a served array.
c--------------------------------------------------------------------------
      
      array_type = array_table(c_array_type,array)

      if (array_type .ne. served_array) then
         print *,'Error: Array type for destroy must be served'
         print *,'op ',(op(i),i=1,loptable_entry)
         print *,'array ' ,array
         print *,'Array table entry: ',(array_table(i,array),
     *     i = 1, larray_table_entry)
         call abort_job()
      endif

      if (me .ne. 0) go to 1000  

c-------------------------------------------------------------------------
c   Build the server message.
c-------------------------------------------------------------------------

      do i = 1, len_sip_server_message
         msg(i) = 0
      enddo

      msg(1) = sip_server_delete_message
      msg(2) = array

c-------------------------------------------------------------------------
c   Send the message to each server.
c-------------------------------------------------------------------------

      do i = 1, nprocs
         company = pst_get_company(i-1)
         if (company .eq. io_company_id) then
            call mpi_send(msg, len_sip_server_message, mpi_integer,
     *               i-1, sip_server_message,
     *               mpi_comm_world, status, ierr)
         endif 
      enddo

 1000 continue

c--------------------------------------------------------------------------
c   Free the persistent blocks of the array.
c--------------------------------------------------------------------------

      call free_persistent_blocks(array, array_table,
     *                    narray_table, index_table, nindex_table,
     *                    block_map_table)

      return
      end
