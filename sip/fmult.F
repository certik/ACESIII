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
      subroutine fmult(array_table, narray_table, 
     *                      index_table,
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
      include 'parallel_info.h'
      include 'epsilon.h'
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
      integer array, index, nindex, ierr
      integer block, seg
      integer find_current_block
      integer blkndx
      integer*8 indblk, get_block_index
      integer stack
      
      integer comm

      integer val1(mx_array_index), val2(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer*8 addr, get_index_from_base

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif
c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      array = op(c_result_array)
      nindex = array_table(c_nindex, array)
      if (nindex .ne. 2) then
         print *,'Error: FMULT argument must be a 2-dimensional array.'
         call abort_job()
      endif

      do i = 1, nindex
         index = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg = index_table(c_current_seg,index)

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

c---------------------------------------------------------------------------
c   Get array data address.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,array) .eq. static_array) then
         addr = address_table(array)
         indblk = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack = array_table(c_array_stack,array)
         indblk = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif 

      if (array_table(c_array_type,array) .eq. static_array) then
         na1 = array_table(c_index_range1, array)
         na2 = array_table(c_index_range2,array)
         ni1 = array_table(c_index_range1+1,array)
         ni2 = array_table(c_index_range2+1,array)
      else
         na1 = val1(1)
         na2 = val2(1)
         ni1 = val1(2)
         ni2 = val2(2)
      endif 

      call eps_mult2(x(indblk), val1(1),val2(1),type(1),
     *                    na1, na2,
     *                    val1(2),val2(2),type(2),ni1,ni2,
     *                    epsilon,epsilonb)    
      
      return
      end

      subroutine eps_mult2(x, i1,i2,itype,ni1,ni2,j1,j2,jtype,
     *                       nj1,nj2,epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer i1,i2,j1,j2
      integer jtype, itype
      integer nj1, nj2, ni1, ni2
      double precision x(ni1:ni2,nj1:nj2)

      double precision epsilon(*), epsilonb(*)

      integer j,i
      double precision eps, epsa, epsi
      double precision val

      if (itype .eq. mobindex) then
         do j = j1, j2
         do i = i1, i2
            x(i,j) = epsilonb(i)*x(i,j)
         enddo
         enddo
      else
         do j = j1, j2
         do i = i1, i2
            x(i,j) = epsilon(i)*x(i,j)
         enddo
         enddo
      endif

      return
      end
