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
      subroutine get_restart_status(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c----------------------------------------------------------------------------
c   Returns 1 or 0 in the array argument, depending on whether the current
c   SIAL program has been restarted or not.  If the SIAL program has been 
c   restarted, a 1 is returned only for the first call to get_restart_status,
c   all other calls to this instruction will return a 0.
c-----------------------------------------------------------------------------
  
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'

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

      integer array, array_type, ind
      double precision val

      array = op(c_result_array)
      if (array .eq. 0) then
         print *,'Error: get_restart_status must be called with a ',
     *           'scalar argument.'
         call abort_job()
      else
         array_type = array_table(c_array_type, array)
         if (array_type .ne. scalar_value) then 
            print *,
     *       'Error: checkpoint instruction requires a scalar ',
     *       'as its argument.'
            call abort_job()
         endif
      endif

      ind =  array_table(c_scalar_index, array)
      if (restart_status) then
         val = 1.0
         restart_status = .false.
      else
         val = 0.0 
      endif

      scalar_table(ind) = val
      return
      end
