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
      subroutine second_moment(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c
c   dipole_moment must be executed with the following syntax:
c   execute dipole_moment x(mu,nu) type ; type = 1, 2, 3, 4, 5, 6
c                                               XX YY ZZ XY XZ YZ
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

      integer i, j, k, type
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
      integer evec_array,nindex_evec,evec_type,sind
      integer a1, a2, b1, b2
      integer m, n, nl2
      integer iatom
      integer*8 iscr 
      integer find_current_block
      integer*8 get_block_data_index
      integer block
      integer igrad, component, intermediate

      double precision dummy, y, sval
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
c   Locate the scalar that tells if derivative of X, Y, or Z
c---------------------------------------------------------------------------
     
      evec_array  = op(c_op1_array)
      evec_type   = array_table(c_array_type, evec_array)
      nindex_evec = array_table(c_nindex, evec_array)

      if (evec_type .ne. scalar_value) return

      if (evec_array .lt. 1 .or. evec_array .gt. narray_table) then
         print *,'Error: Invalid scalar in dipole_moment, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()          
      endif
         
      sind =  array_table(c_scalar_index, evec_array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
        print *,'Scalar table index out of range in dipole_moment, ',
     *           'line ',current_line
        print *,'Index for array ',evec_array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      sval = scalar_table(sind)

      type = 0
      if ( (abs (sval - 1.0)) .lt. 0.001 ) type = 1
      if ( (abs (sval - 2.0)) .lt. 0.001 ) type = 2
      if ( (abs (sval - 3.0)) .lt. 0.001 ) type = 3
      if ( (abs (sval - 4.0)) .lt. 0.001 ) type = 4
      if ( (abs (sval - 5.0)) .lt. 0.001 ) type = 5
      if ( (abs (sval - 6.0)) .lt. 0.001 ) type = 6


c---------------------------------------------------------------------------
c   Compute data block of integrals.
c---------------------------------------------------------------------------

      if (intpkg .eq.flocke_package) then
         call comp_second_moment(type,x(indblk), a1, a2, b1, b2,nsend)
      else
        write(6,*) ' MUST USE ERD IN FOCK_DER!!! ' 
      endif 

      return
      end
