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
      subroutine set_index(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c---------------------------------------------------------------------------
c   Sets the indices of a 4-d static array in common block values.
c   The four indeces will be used as flags indicating which index to 
c   truncate. These indices are stored in the SINDEX common block. 
c
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      common /SINDEX/index1, index2, index3, index4 
      integer index1, index2, index3, index4  
      double precision flags_value

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

      integer ierr, array, array_type, ind, nind
      integer i

      array = op(c_result_array)
      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: Invalid array in set_flags, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif

      nind = array_table(c_nindex, array)
      if (nind .ne. 4) then
         print *,'Error: set_index requires a 4-index array.'
         call abort_job()
      endif

c-----------------------------------------------------------------------
c   Index1-4 are determined from the c_current_seg
c   field of the 1st, 2nd, 3rd and 4th index of the array.
c------------------------------------------------------------------------

      ind    = array_table(c_index_array1,array)
      index1 = index_table(c_current_seg, ind)

      ind    = array_table(c_index_array1+1,array)
      index2 = index_table(c_current_seg, ind)

      ind    = array_table(c_index_array1+2,array)
      index3 = index_table(c_current_seg, ind)

      ind    = array_table(c_index_array1+3,array)
      index4 = index_table(c_current_seg, ind)

      if (dbg) 
     *  write(6,*) 'Task ',me,' RINDEX :', index1, index2, 
     *              index3, index4 

      return
      end
