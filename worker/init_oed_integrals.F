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
      subroutine init_oed_integrals(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table, 
     *                      address_table, op, type)
c---------------------------------------------------------------------------
c   Generates one-electron integrals of a particular type in the memory
c   allocated for the array in the array of "op".
c
c   type = 1  --> NAI integrals
c          2  --> KIN integrals
c          3  --> OVL integrals
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
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
      integer*8 address_table(narray_table)
      integer type

      integer array
      integer*8 addr, ixarray, get_index_from_base
      double precision xxx(1)
#ifdef ALTIX
      pointer (dptr, xxx)
#else
      common xxx
#endif

c-------------------------------------------------------------------------
c   Locate array address in array_table.
c-------------------------------------------------------------------------

#ifdef ALTIX
      dptr = dshptr
#endif

      array    = op(c_result_array)
      addr     = address_table(array)
      ixarray  = get_index_from_base(addr, xxx, 2)

      call compute_oed_integrals(type, xxx(ixarray))

      return
      end
