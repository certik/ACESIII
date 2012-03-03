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
      subroutine smon_on(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
      implicit none
      include 'interpreter.h'
      include 'server_monitor.h'
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
      double precision scalar_table(nscalar_table) 
      integer*8 address_table(narray_table)

      logical is_file_open
      character*128 fn

      server_monitor_on = .true.

      server_logunit = 89
      inquire (unit=server_logunit, opened = is_file_open)
      if (is_file_open) close (server_logunit)

      call generate_scratch_filename(fn, 'smon')
      open (file = fn, unit = server_logunit)
      write (server_logunit,*) 'SERVER MONITOR at line ',current_line

      return
      end
