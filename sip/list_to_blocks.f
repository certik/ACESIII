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
      subroutine list_to_blocks(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table,  op)
c--------------------------------------------------------------------------
c   Reads block data from a list file and distributes them to their owners.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h'
      include 'saved_data.h'

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

      integer i, j, k, n
      integer array, index(mx_array_index), nindex, ierr
      integer block, seg(mx_array_index)
      integer blkndx, get_block_number
      integer saveseg(mx_array_index)
      integer*8 indblk, get_block_data_index
      integer stack
      integer f_form_msg_tag
      integer datatag
      integer type, kreq, ireq, nreq
      
      integer mx_blk
      parameter (mx_blk = 25000)
      integer request(mx_blk)
      integer drequest(mx_blk)
      integer desc_msg(2,mx_blk)
      integer statuses(MPI_STATUS_SIZE,mx_blk)

      integer my_comm_rank, comm
      integer home, iblk
      integer val1(mx_array_index), val2(mx_array_index)
    
      integer nb, proc, map, company_comm
      integer tag
      integer pst_get_company_comm, pst_get_company_rank
      integer status(MPI_STATUS_SIZE)

      if (first_list_to_blocks) then
         first_list_to_blocks = .false.
         narray_list          = 0
      endif

      narray_list = narray_list + 1
      if (narray_list .gt. max_array_list) then
         print *,'Error: LIST_TO_BLOCKS arrays exceeds maximum of ',
     *       max_array_list
         call abort_job()
      endif

      array_list(narray_list) = op(c_result_array)
      return
      end
