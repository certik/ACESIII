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
      subroutine compute_diis(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Establishes the lower triangular array elements for a DIIS calculation.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
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

      integer i, handle, ind

      integer diis_counter, bhandle
      double precision b_array
      common /diis/b_array(10,10),
     *               diis_counter, bhandle(10)

      double precision c1, c2, c3, c4, c5, c6, c7, c8, c9, c10

c-------------------------------------------------------------------------
c   Put the array index in the appropriate slot.
c-------------------------------------------------------------------------

      if (diis_counter .ne. 55) then
         print *,'Error: COMPUTE_DIIS at line ',current_line
         print *,'DIIS_SETUP has only been called ',diis_counter,
     *      ' times, should be 55.'
         call abort_job()
      endif

      call form_R(b_array, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 )

c-------------------------------------------------------------------------
c   Store the "c" values back in the diagonal elements of "b".
c-------------------------------------------------------------------------

      ind = array_table(c_scalar_index,bhandle(1))
      scalar_table(ind) = c1

      ind = array_table(c_scalar_index,bhandle(2))
      scalar_table(ind) = c2

      ind = array_table(c_scalar_index,bhandle(3))
      scalar_table(ind) = c3

      ind = array_table(c_scalar_index,bhandle(4))
      scalar_table(ind) = c4

      ind = array_table(c_scalar_index,bhandle(5))
      scalar_table(ind) = c5

      ind = array_table(c_scalar_index,bhandle(6))
      scalar_table(ind) = c6

      ind = array_table(c_scalar_index,bhandle(7))
      scalar_table(ind) = c7

      ind = array_table(c_scalar_index,bhandle(8))
      scalar_table(ind) = c8

      ind = array_table(c_scalar_index,bhandle(9))
      scalar_table(ind) = c9

      ind = array_table(c_scalar_index,bhandle(10))
      scalar_table(ind) = c10

      return
      end

