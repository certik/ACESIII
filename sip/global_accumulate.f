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
      subroutine global_accumulate(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, address_table,
     *                      comm, op, comm_timer, instruction_timer)
c---------------------------------------------------------------------------
c   If the current block of the result array is mapped locally,   
c   global_accumulate sums the operand array into the local block.
c
c   If the current block is mapped to another processor, we perform an 
c   accumulate operation.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      include 'trace.h'

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      integer*8 address_table(narray_table)
      double precision scalar_table(nscalar_table)

      integer result_array, op1_array, dest
      integer result_block, op1_block
      integer find_block_destination
      integer ierr
      integer i
      integer flopcount, send_count
      double precision send_time
      integer find_current_block
      integer blkndx
      integer find_current_block_map
      integer comm
      integer local_block, result_blkndx
      integer op1_type
      integer comm_timer, instruction_timer

      result_array = op(c_result_array)
      op1_array    = op(c_op1_array)
      op1_type     = array_table(c_array_type,op1_array)
      if (op1_type .eq. static_array) then
         print *,'Error: Attempt to use a static block in an ',
     *           'accumulate instruction.'
         print *,'Local array is array ',op1_array
         call abort_job() 
      endif

c-------------------------------------------------------------------------
c   Locate the current block of the operand array.
c-------------------------------------------------------------------------

      op1_block = find_current_block(op1_array, 
     *                             array_table(1,op1_array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)
      if (op1_block .le. 0) then
         print *,'Error in global_accumulate on task ',me,
     *          ' Cannot find block for op1_array = ',op1_array
         print *,'Line ',current_line
         print *,'current segs = ',(index_table(c_current_seg,i),
     *        i=1,nindex_table)
         call dump_block_ids()
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Locate the current block of the result array.
c-------------------------------------------------------------------------

      result_block = find_current_block_map(result_array, 
     *                                  array_table(1,result_array),
     *                                  index_table, nindex_table,
     *                                  block_map_table,
     *                                  nblock_map_table)

c--------------------------------------------------------------------------
c   Determine location of the mapping of the result_array's current block.
c--------------------------------------------------------------------------

      dest = find_block_destination(result_array, array_table,
     *               narray_table, index_table, nindex_table,
     *               block_map_table, nblock_map_table)
 
      if (dest .eq. my_company_rank) then

c--------------------------------------------------------------------------
c   Result block must exist already.
c--------------------------------------------------------------------------

         if (result_block .le. 0) then
            print *,'Error in global_accumulate on task ',me,
     *          ' Cannot find block for result_array = ',result_array
            print *,'current segs = ',(index_table(c_current_seg,i),
     *        i=1,nindex_table) 
            print *,'dest, me, my_company_rank = ',
     *            dest, me, my_company_rank
            call dump_block_ids()
            call abort_job()
         endif

         local_block = find_current_block(result_array,
     *                             array_table(1,result_array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, result_blkndx)
      else
         result_blkndx = 0
      endif

c---------------------------------------------------------------------------
c   Accumulate into the block at the destination processor.
c---------------------------------------------------------------------------

      call accumulate_block(op1_array, op1_block, blkndx,
     *                 result_array, result_block, result_blkndx,
     *                 dest, array_table,
     *                 narray_table, index_table, nindex_table,
     *                 segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table,
     *                 scalar_table, nscalar_table, address_table,
     *                 send_count, send_time, comm, my_company_rank,
     *                 op(c_opcode), comm_timer, instruction_timer) 
      return
      end
