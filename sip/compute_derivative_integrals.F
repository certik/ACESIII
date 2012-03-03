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
      subroutine compute_derivative_integrals(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Multiplies the each sample of the block given by the array argument by 
c   diagonal eigenvalue matrix.  The eigenvalues matrix is not stored in
c   a 2-dimensional form, since it can be represented by a vector of
c   the diagonal elements.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
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

      integer i, j, k
      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer array, index, nindex, ierr
      integer blk, blkndx, seg
      integer nsend
      integer*8 indblk(12), get_block_data_index
      integer der_flags(12)
      integer stack
      
      integer me, my_comm_rank, comm

      integer val1(mx_array_index,12), val2(mx_array_index,12)
      double precision x(1)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (iptr, ibuf)
      pointer (xptr, x)
#else
      equivalence (ibuf, dbuf)
      common dbuf
#endif

      integer flag, iarray, handle

      integer dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,
     *        derint_counter
      integer derint_array(12)
      common /derint/dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,
     *               derint_counter
      equivalence (derint_array(1), dx1)

      if (derint_counter .ne. 12) then
         print *,'Error: compute_derivative_integrals mapping ',
     *     'is not set properly. Number of arrays is ',derint_counter
         call abort_job()
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Create a block for each of the arrays.
c---------------------------------------------------------------------------

      do iarray = 1, 12
         handle = derint_array(iarray)
         call create_current_block(handle,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op,
     *                 .true., .false., blk, ierr)
         blkndx = ierr

         call get_block_computed_flag(handle, blk, blkndx, flag)
         if (flag .eq. 0) then
            call set_opblock(handle, blk, blkndx, op)
            call set_block_computed_flag(handle, blk, blkndx, 1)
         endif
                                                                                
         stack = array_table(c_array_stack,handle)
         indblk(iarray) = get_block_data_index(handle, blk, stack,
     *                                         blkndx, x)

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
         nindex = array_table(c_nindex, handle)

         do i = 1, nindex
            index = array_table(c_index_array1+i-1,handle)
            seg = index_table(c_current_seg,index)

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

            call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i,iarray), 
     *                             val2(i,iarray))
         enddo   ! index loop
      enddo      ! array loop

c--------------------------------------------------------------------------
c   Block ranges must match to be valid.
c--------------------------------------------------------------------------

      a1 = val1(1,1)
      a2 = val2(1,1)
      b1 = val1(2,1)
      b2 = val2(2,1)
      c1 = val1(3,1)
      c2 = val2(3,1)
      d1 = val1(4,1)
      d2 = val2(4,1)

      do i = 2, 12
         if (val1(1,i) .ne. a1 .or. val2(1,i) .ne. a2 .or.
     *       val1(2,i) .ne. b1 .or. val2(2,i) .ne. b2 .or.
     *       val1(3,i) .ne. c1 .or. val2(3,i) .ne. c2 .or.
     *       val1(4,i) .ne. d1 .or. val2(4,i) .ne. d2) then
            print *,'Error: Indices for all der_int_setup arrays ',
     *              'must match.'
            print *,'Array ',i,' does not match array 1.'
            call abort_job()
         endif
      enddo

c-----------------------------------------------------------------------------
c   Calculate the integral blocks.
c-----------------------------------------------------------------------------

      if (intpkg .eq. flocke_package) then
         do i = 1, 12
            do j = 1, 12
               der_flags(j) = 0
            enddo

            der_flags(i) = 1
            call compute_erd_2e_derv_integrals(a1,a2,b1,b2,c1,c2,d1,d2,
     *                 der_flags,
     *                 dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *                 dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *                 ibuf(iccbeg), ibuf(iccend),
     *                 x(indblk(i)), nsend)
         enddo
      else
         print *,'Error: Invalid integral package for',
     *            ' derivative integrals: ',intpkg
         call abort_job() 
      endif

      return
      end
