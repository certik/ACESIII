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
      subroutine accumulate_block(array, array_block, iblk2,
     *                 tarray, tarray_block, iblk, tproc, array_table,
     *                 narray_table, index_table, nindex_table,
     *                 segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table,
     *                 scalar_table, nscalar_table, address_table,
     *                 send_count, send_time, comm, my_comm_rank,
     *                 opcode, comm_timer, instruction_timer)
c--------------------------------------------------------------------------
c   Sends a block of a distributed array to its destination in rsponse to 
c   an accumulate command.  
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blkmgr.h'

      integer array, array_block
      integer tarray, tarray_block
      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer*8 address_table(narray_table)
      integer send_count, flopcount
      integer comm
      integer my_comm_rank
      integer opcode

      integer result_type, result_block, tproc, flag
      integer array_type
      integer find_block_destination
      integer iblk, iblk2, get_block_number
      integer nwblock, request
      integer get_block_request
      integer comm_timer, instruction_timer
      integer tblkndx

      double precision send_time
      integer nscalar_table
      double precision scalar_table(nscalar_table)

      result_type = array_table(c_array_type,tarray)
      array_type = array_table(c_array_type, array)
      if (result_type .ne. distributed_array) then
         print *,'Error: accumulate_block called, ',
     *             ' but the array is not a distributed array.'
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Set the "put" flag, indicating that the barrier should destroy
c   all local copies of the array's blocks.
c--------------------------------------------------------------------------

      array_table(c_put_flag,tarray) = 1

      call mutex_block(iblk2)   ! lock the "send" block

c---------------------------------------------------------------------------
c   Check for a previous communication request.  It must be satisfied before
c   a new one can be started.
c---------------------------------------------------------------------------

      request = get_block_request(array, array_block, iblk2)
      if (request .ne. mpi_request_null) 
     *   call wait_on_block(array, array_block, iblk2, 
     *             array_type, request,
     *             instruction_timer, comm_timer)

      if (tproc .eq. my_comm_rank) then
         call mutex_block(iblk)

         if (opcode .eq. put_op) then

c--------------------------------------------------------------------------
c   Sum the block locally.
c--------------------------------------------------------------------------

            call sum_blocks(tarray, tarray_block, iblk, 
     *                      array, array_block, iblk2,
     *                      tarray, tarray_block, iblk, sum_op,
     *                      array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
         else if (opcode .eq. put_replace_op) then

c----------------------------------------------------------------------------
c   Copy the block locally.
c----------------------------------------------------------------------------

            call assign_block(array, array_block, iblk2,
     *                      tarray, tarray_block, iblk,
     *                      array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
         else
            print *,'Error: accumulate_block must be called with ',
     *              'put_op or put_replace_op.'
            print *,'Invalid opcode was ',opcode
            call abort_job()
         endif

         call release_mutex_block(iblk)
      else

c--------------------------------------------------------------------------
c   Send the data to a remote process.
c--------------------------------------------------------------------------

         call get_actual_blocksize(array, array_block, iblk2,
     *                  array_table, narray_table, 
     *                  index_table, nindex_table, 
     *                  segment_table, nsegment_table, nwblock)

         tblkndx = -1

         call send_block(array, array_block, iblk2, tproc, tarray, 
     *                tarray_block, tblkndx, nwblock, comm, opcode,
     *                request)

         call set_block_request(array, array_block, iblk2, request)
         if (request .ne. mpi_request_null) 
     *      call blkmgr_insert_block_in_list(comm_list_head, 
     *                     comm_list_tail, iblk2,
     *                     c_comm_list_ptr, .true.)
       endif

      call release_mutex_block(iblk2)  
      return
      end

