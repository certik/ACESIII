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
      subroutine collective_sum(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      comm, op)
c---------------------------------------------------------------------------
c   Handler for the collective_sum_op operation.
c   This routine performs a global accumulate of a scalar on all processors,
c   (i. e. mpi_all_reduce).
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table, comm
      double precision scalar_table(nscalar_table)

      integer i, result_array, op1_array
      integer iop1, ires
      integer nprocs
      integer ierr

      result_array = op(c_result_array)
      op1_array    = op(c_op1_array)

c-------------------------------------------------------------------------
c   Check for scalar operands.
c-------------------------------------------------------------------------

      if (array_table(c_array_type,result_array) .ne. 
     *                              scalar_value .or.
     *    array_table(c_array_type,op1_array) .ne. scalar_value) then
         print *,'Error in collective_sum: Arrays are not scalars.'
         print *,'Op1 = ',op1_array,' Result = ',result_array
         print *,'Operation: ',(op(i),i=1,loptable_entry)
         call abort_job()
      endif 

c-------------------------------------------------------------------------
c   Look up the scalar values in the scalar_table.
c-------------------------------------------------------------------------

      iop1 = array_table(c_scalar_index,op1_array)
      ires = array_table(c_scalar_index,result_array)

c-------------------------------------------------------------------------
c   Accumulate the sum into the result value.
c-------------------------------------------------------------------------

      call mpi_comm_size(comm, nprocs, ierr)
      if (nprocs .eq. 1) then
         scalar_table(ires) = scalar_table(ires) + scalar_table(iop1)
      else   
         call mpi_allreduce(scalar_table(iop1), scalar_table(ires), 
     *              1, mpi_double_precision, mpi_sum, comm, ierr)
      endif
      return
      end
