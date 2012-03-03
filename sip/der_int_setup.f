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
      subroutine der_int_setup(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Establishes the array indices for a derivative integral calculation.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
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

      integer i, j, k
      integer flag, handle, index, nindex, ierr
      integer blk, blkndx

      integer dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,
     *        derint_counter
      integer derint_array(12)
      common /derint/dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,
     *               derint_counter
      equivalence (derint_array(1), dx1)

      if (first_der_int_setup) then
         derint_counter = 0
         first_der_int_setup = .false.
      endif
      
c-------------------------------------------------------------------------
c   Put the array index in the appropriate slot.
c-------------------------------------------------------------------------

      derint_counter = derint_counter + 1
      if (derint_counter .gt. 12) derint_counter = 1

      handle = op(c_result_array)
      derint_array(derint_counter) = handle

c------------------------------------------------------------------------
c   Make sure the requested block exists.  create_current_block will
c   create the block if it does not exist, and simply return if the
c   block is already present.
c------------------------------------------------------------------------
                                                                                
      call create_current_block(handle,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., .false., blk, ierr)
      blkndx = ierr
                                                                                
      call get_block_computed_flag(handle, blk, blkndx, flag)
      if (flag .eq. 0) then
         call set_opblock(handle, blk, blkndx, op)
         call set_block_computed_flag(handle, blk, blkndx, 1)
      endif

      return
      end

