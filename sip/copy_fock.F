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
      subroutine copy_fock(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Copies one array into a Fock array and copies the diagonal elements
c   into the corresponding eigenvalue array.
c
c   execute Array1 fock_a 
c         OR
c   execute Array1 fock_b
c 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
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

      integer*8 iarray, ifock, get_index_from_base
      integer i, j, n
      integer array, fock_array
      integer nindex
      integer ind(mx_array_index)

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

c----------------------------------------------------------------------------
c   Locate the data for the input matrix.
c---------------------------------------------------------------------------

#ifdef ALTIX
      dptr = dshptr
#endif

      array = op(c_result_array)
      fock_array = op(c_op1_array)
      if (array .eq. 0 .or. fock_array .eq. 0) then
         print *,'Error: COPY_FOCK routine requires 2 array arguments.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      if (fock_array .ne. 9 .and. 
     *    fock_array .ne. 10) then
         print *,'Error: COPY_FOCK requires a Fock array as arg. 2'
         call abort_job()
      endif

      nindex = array_table(c_nindex, array)

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,fock_array)
      enddo

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray = get_index_from_base(address_table(array), x, 2)
      ifock  = get_index_from_base(address_table(fock_array), x, 2)

      n = index_table(c_index_size, ind(1))  ! pick up length of index
      
      do i = 1, n*n
         x(ifock+i-1) = x(iarray+i-1)
      enddo

c-------------------------------------------------------------------------
c   Copy diagonal elements into epsilon or epsilonb.
c-------------------------------------------------------------------------

      if (index_table(c_index_type,ind(1)) .eq. moaindex) then
         call copy_diag(x(ifock), n, epsilon)
      else
         call copy_diag(x(ifock), n, epsilonb)
      endif

      return
      end

      subroutine copy_diag(x, n, y)
      implicit none
      integer i, n
      double precision x(n,n), y(n)

      do i = 1, n
         y(i) = x(i,i)
      enddo
      return
      end
