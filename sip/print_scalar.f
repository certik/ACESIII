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
      subroutine print_scalar(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c---------------------------------------------------------------------------
c   Prints the value of a scalar variable.  The scalar to be printed
c   is defined in the c_result_array field of the op argument.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
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

      integer ierr, array, array_type, ind
      integer i

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      if (array_type .ne. scalar_value) return

      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: Invalid array in print_scalar, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif

      ind =  array_table(c_scalar_index, array)
      if (ind .lt. 1 .or. ind .gt. nscalar_table) then
         print *,'Scalar table index out of range in print_scalar, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      if (dbg) then
         print *,'Task ',me,' Scalar #',array,': value = ',
     *      scalar_table(ind),' at line number ',current_line
         call prt_time('Worker time')
      else
         if (me .eq. 0) 
     *      print *,'Task ',me,' Scalar #',array,': value = ',
     *      scalar_table(ind),' at line number ',current_line
      endif
      call c_flush_stdout()
      return
      end
