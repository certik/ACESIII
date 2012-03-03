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
      subroutine scf_rhf(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c-----------------------------------------------------------------------------
c   This routine copies the current values of the CA pre-defined array
c   into the CB array, and the FOCKA int FOCKB.  This effectively duplicates 
c   the alpha-spin
c   eigenvectors and eigenvalues as beta-spin.  It is used only in the
c   RHF versions of SCF codes, enabling UHF codes to be used from this
c   point in the processing (i. e. tran, ccsd, lambda, ...).
c-----------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'trace.h'
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

      integer*8 get_index_from_base
      integer*8 ixca, ixcb, ixfocka, ixfockb
      integer j, k, n, index, bseg, eseg, val1, val2
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)

      dptr = dshptr
#else
      common x
#endif

c----------------------------------------------------------------------------
c   Find the address of each array. (CA is array 2, CB is array 3, FOCKA is
c   array 9, FOCKB is array 10).
c----------------------------------------------------------------------------

      ixca = get_index_from_base(address_table(2), x, 2)
      ixcb = get_index_from_base(address_table(3), x, 2)
      ixfocka = get_index_from_base(address_table(9), x, 2)
      ixfockb = get_index_from_base(address_table(10), x, 2)

c----------------------------------------------------------------------------
c   Determine size of the arrays.
c----------------------------------------------------------------------------

      n = 1
      do j = 1, array_table(c_nindex,2)
         index = array_table(c_index_array1+j-1,2)
         bseg = index_table(c_bseg, index)
         eseg = index_table(c_eseg, index)
         call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
         array_table(c_index_range1+j-1,2) = val1
         call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
         array_table(c_index_range2+j-1,2) = val2
         k = array_table(c_index_range2+j-1,2) -
     *       array_table(c_index_range1+j-1,2) + 1
         n = n * k
      enddo

c---------------------------------------------------------------------------
c   Copy the data.
c---------------------------------------------------------------------------

      do j = 1, n
         x(ixcb+j-1) = x(ixca+j-1)
      enddo

c----------------------------------------------------------------------------
c   Copy the FOCKA into FOCKB.
c----------------------------------------------------------------------------

      do j = 1, n
         x(ixfockb+j-1) = x(ixfocka+j-1)
      enddo
      return
      end
