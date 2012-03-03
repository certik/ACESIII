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
      subroutine broadcast_array(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c---------------------------------------------------------------------------
c   Broadcasts the static array in the "result" field of the instruction
c   to each processor in the company.  The copy of the array on proc 0
c   is the data that is broadcast.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

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

      integer ierr, array, array_type, ind
      integer i

      double precision x(1)
#ifdef ALTIX
      pointer (ptr, x)
#else
      common x
#endif

      integer*8 daddr, indx
      integer*8 get_index_from_base
      integer n, nindex, val1, val2
      integer company_comm, pst_get_company_comm

#ifdef ALTIX
      ptr = dshptr
#endif
      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      if (array_type .ne. static_array) then
         print *,'Error in broadcast_array: Arg must be a static array'
         call abort_job()
      endif

      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: Invalid array in get_my_rank, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Locate the array's memory and determine the size of the array.
c---------------------------------------------------------------------------

      daddr = address_table(array)
      indx = get_index_from_base(daddr, x, 2)
      nindex = array_table(c_nindex,array)

      n = 1
      do i = 1, nindex
         val1 = array_table(c_index_range1+i-1,array)
         val2 = array_table(c_index_range2+i-1,array)
         n = n * (val2-val1+1)
      enddo

c----------------------------------------------------------------------------
c   Broadcast the data.
c----------------------------------------------------------------------------

      company_comm = pst_get_company_comm(me)
      call mpi_bcast(x(indx), n, MPI_DOUBLE_PRECISION, 0,
     *               company_comm, ierr)
      return
      end
