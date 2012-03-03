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
      subroutine compute_abac_batch(array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c-------------------------------------------------------------------------- 
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, array 
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table), get_index_from_base
      integer*8 get_block_data_index
      integer find_current_block, block 

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 indblk, get_block_index
      integer*8 integral_scr
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack

      double precision x(1)
      double precision dbuf(1), temp 
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (xptr, x)
      pointer (iptr, ibuf)
#else
      common dbuf, x
      equivalence (ibuf(1), dbuf(1))
#endif

      array = op(c_result_array)
      direct_flag = .false. 

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c------------------------------------------------------------------------
c   Make sure the requested block exists.  create_current_block will
c   create the block if it does not exist, and simply return if the 
c   block is already present.
c------------------------------------------------------------------------

      call create_current_block(array,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., direct_flag, blk, ierr)
      blkndx = ierr
      call get_block_computed_flag(array, blk, blkndx, flag)
      if (flag .eq. 0) then
         call set_opblock(array, blk, blkndx, op)
         call set_block_computed_flag(array, blk, blkndx, 1)
      endif

      stack = array_table(c_array_stack,array)
      indblk = get_block_index(array, blk, stack,
     *                         blkndx, x, .true.)

c--------------------------------------------------------------------------
c   Get the segment ranges for each index of the integral array.
c--------------------------------------------------------------------------

      nind = array_table(c_nindex,array)
      if (nind .ne. 4) then
         print *,'Error: local_integral_block was called with array ',
     *        array,' which has ',nind,' indices.'
         print *,'       This routine requires 4 indices.'
         call abort_job()
      endif

      do i = 1, nind
         ind = array_table(c_index_array1+i-1,array)
         seg = index_table(c_current_seg,ind)
         call get_index_segment(ind, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

c--------------------------------------------------------------------------
c   Form the address of the current data block.
c--------------------------------------------------------------------------

      block = find_current_block(array, array_table(1,array),
     *                          index_table, nindex_table,
     *                          segment_table, nsegment_table,
     *                          block_map_table, blkndx)

      stack = array_table(c_array_stack, array)
      indblk = get_block_data_index(array, block, stack,
     *                              blkndx, x)

c---------------------------------------------------------------------------
c   Compute the block of integrals.
c---------------------------------------------------------------------------

         if (intpkg .eq. flocke_package) then
            integral_scr = dpkgscr
            maxblk  = zmax
         else if (intpkg .eq. gamess_package .or.
     *            intpkg .eq. gamess_derivative_package) then
            integral_scr = i0wrk
            maxblk = 15**4    ! max. angular momentum for GAMESS is 15.
         endif

      if (.not. dryrun) 
     *   call compute_abac_integrals(val1(1),val2(1), val1(2),val2(2),
     *                       val1(3),val2(3),val1(4),val2(4),
     *                 dbuf(integral_scr),maxblk, ibuf(ipkgscr),
     *                 dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha), 
     *                 ibuf(iccbeg), ibuf(iccend), 
     *                 x(indblk), nsend)

      return
      end
