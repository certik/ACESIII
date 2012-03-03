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
      subroutine array_copy(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Broadcast a array_copy message to all servers in the IOCOMPANY.
c   Format of command is :
c
c   execute array_copy source target
c
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
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, ii, source, target, company, comm, ierr
      integer source_type, target_type
      integer temp
      integer msg(len_sip_server_message)
      integer status(MPI_STATUS_SIZE) 
      integer pst_get_company
      integer pst_get_company_comm

      target = op(c_op1_array)
      source = op(c_result_array)
      comm = pst_get_company_comm(me)

c--------------------------------------------------------------------------
c   Pick up the array indices.  Both array types must be served_array.
c   The instruction requires 2 arrays: the source array and the target
c   array.
c--------------------------------------------------------------------------
      
      source_type = array_table(c_array_type,source)
      target_type = array_table(c_array_type,target)

      if (source_type .ne. target_type) then
         print *,'Error: Array types in array_copy do not match.'
         print *,'Array ',source,' array_type is ',source_type,
     *     ' target array ',target,' array_type is ',
     *     target_type
         call abort_job()
      endif
 
      if (source_type .ne. served_array) then
         print *,'Error: Array type for array_copy must be served'
         call abort_job()
      endif

      if (me .ne. 0) go to 1000  

c-------------------------------------------------------------------------
c   Build the server message.
c-------------------------------------------------------------------------

      do i = 1, len_sip_server_message
         msg(i) = 0
      enddo

      msg(1) = sip_server_copy_message
      msg(2) = target
      msg(3) = source

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
c   Free the persistent blocks of both the source and target arrays.
c--------------------------------------------------------------------------

      call free_persistent_blocks(source, array_table,
     *                    narray_table, index_table, nindex_table,
     *                    block_map_table)

      call free_persistent_blocks(target, array_table,
     *                    narray_table, index_table, nindex_table,
     *                    block_map_table)

c-------------------------------------------------------------------------
c   Now swap the array table entries.  Since one of these entries is a
c   pointer into the block_map_table, this in effect causes all the 
c   blocks of each array to be remapped.
c-------------------------------------------------------------------------

      do i = 1, larray_table_entry
         temp = array_table(i,source)
         array_table(i,source) = array_table(i,target)
         array_table(i,target) = temp
      enddo

      return
      end
