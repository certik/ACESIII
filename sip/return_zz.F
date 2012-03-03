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
      subroutine return_zz(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c
c   return_zz must be executed with the following syntax:
c   execute return_zz x(mu,nu)
c
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
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

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 indblk, get_block_index
      integer*8 integral_scr
      integer val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
      integer array, ind(mx_array_index), seg(mx_array_index)
      integer a1, a2, b1, b2
      integer m, n, nl2
      integer iatom
      integer*8 iscr 
      integer find_current_block
      integer*8 get_block_data_index
      integer block
      integer igrad, component, intermediate

      double precision dummy, y
      integer flopcount
      integer*8 arg64(10)
      logical*8 l8false

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter, dcoord

      nl2 = (nbasis * nbasis + nbasis)/2
      handle = op(c_result_array)
      direct_flag = .false.

#ifdef ALTIX
      dptr = dshptr
#endif

c------------------------------------------------------------------------
c   Make sure the requested block exists.  create_current_block will
c   create the block if it does not exist, and simply return if the
c   block is already present.
c------------------------------------------------------------------------

      call create_current_block(handle,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., direct_flag, blk, ierr)
      blkndx = ierr
      call get_block_computed_flag(handle, blk, blkndx, flag)
      if (flag .eq. 0) then
         call set_opblock(handle, blk, blkndx, op)
         call set_block_computed_flag(handle, blk, blkndx, 1)
      endif

      stack = array_table(c_array_stack,handle)
      indblk = get_block_index(handle, blk, stack,
     *                         blkndx, x, .true.)

c--------------------------------------------------------------------------
c   Find the ranges of the AO indices of the input block.
c--------------------------------------------------------------------------

         array = op(c_result_array)

         nind  = array_table(c_nindex,array)
         if (nind .ne. 2) then
            print *,'Error: Fock_der requires a 2-index argument array.'
            print *,'Array ',array,' is defined with ',nind,' indices.'
            call abort_job()
         endif       
   
         do i = 1, nind
            ind(i) = array_table(c_index_array1+i-1,array)
            seg(i) = index_table(c_current_seg,ind(i))
         enddo

         call get_index_segment(ind(1), seg(1), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, a1, a2)
         call get_index_segment(ind(2), seg(2), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, b1, b2)

c--------------------------------------------------------------------------
c   Form the address of the current data block.
c--------------------------------------------------------------------------

         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack = array_table(c_array_stack, array)
         indblk = get_block_data_index(array, block, stack,
     *                                        blkndx, x)

c---------------------------------------------------------------------------
c   Compute data block of integrals.
c---------------------------------------------------------------------------

      if (intpkg .eq.flocke_package) then
         call comp_return_zz(x(indblk), a1, a2, b1, b2,nsend)
      else
        write(6,*) ' MUST USE ERD IN FOCK_DER!!! ' 
      endif 

      return
      end
