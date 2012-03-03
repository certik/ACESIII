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
      subroutine remove_ds(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)

c--------------------------------------------------------------------------
c
c   
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
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
      integer type(mx_array_index) 
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer*8 iarray, ievec, get_index_from_base
      integer i, j, n1, n2, n3, n4  
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer ind(mx_array_index)
      integer seg, index(mx_array_index), val1(mx_array_index), 
     *        val2(mx_array_index), junk

      integer block, blkndx
      integer find_current_block
      integer*8 get_block_index
      integer stack

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
#ifdef ALTIX
      dptr = dshptr
#endif

c----------------------------------------------------------------------------
c   Locate the data for the input matrix.
c---------------------------------------------------------------------------

      array      = op(c_result_array)
      evec_array = op(c_op1_array)

      if (array .eq. 0 .or. evec_array .eq. 0) then
         print *,'Error: remove_single routine requires 2 
     *            array arguments.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      array_type = array_table(c_array_type, array)
      evec_type  = array_table(c_array_type, evec_array)

      nindex      = array_table(c_nindex, array)
      nindex_evec = array_table(c_nindex, evec_array)

      if (nindex .ne. 2) then 
         write(6,*) ' Remove_ds must use arrays with two indices' 
         write(6,*) ' Number given =', nindex, nindex_evec 
         call abort_job() 
      endif 

      if (nindex .ne. nindex_evec) then
         print *,'Error: Both arrays in return_diagonal must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex,' indices.'
         print *,'Second array has ',nindex_evec,' indices.'
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      block  = find_current_block(array, array_table(1,array),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,array) 
      iarray = get_block_index(array, block, stack, blkndx, x, .true.)  

      block  = find_current_block(evec_array, array_table(1,evec_array),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,evec_array) 
      ievec = get_block_index(evec_array, block, stack, blkndx,x,.true.)  

c     iarray = get_index_from_base(address_table(array), x, 2)
c     ievec  = get_index_from_base(address_table(evec_array), x, 2)

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------

      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1,array)
         type(i)  = index_table(c_index_type, index(i))
         seg      = index_table(c_current_seg,index(i))

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

         call get_index_segment(index(i), seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo


      call retnods(x(iarray),x(ievec),
     *              val1(1),val2(1),val1(2),val2(2),index,type)

      return 
      end 

      subroutine retnods(array1,array2,a1,a2,b1,b2,index,type) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'  

      integer i, a, b, c, d, a1, a2, b1, b2, c1, c2, d1, d2  
      double precision array1(a1:a2,b1:b2) 
      double precision array2(a1:a2,b1:b2) 

      integer index(mx_array_index), type(mx_array_index) 
      integer nzero, zero(100), index1, index2, index3, index4 
      integer nszero, ndzero, szero(100), dzero(100) 
      common /SINDEX/index1, index2, index3, index4 

c     Initialize output array to input array. 
c     --------------------------------------- 

      do a = a1, a2 
      do b = b1, b2 

         array2(a,b) = array1(a,b) 

      enddo ! b 
      enddo ! a 
c
c     Determine range of singly occupied indices to be removed. 
c     --------------------------------------------------------- 

      b      = 0 
      nszero = 0 
      if (nalpha_occupied .gt. nbeta_occupied) then 
         nszero = nalpha_occupied - nbeta_occupied 
         do a = nalpha_occupied, nbeta_occupied, -1 
            b = b + 1 
            szero(b) = a 
         enddo 
      endif 
      if (nbeta_occupied .gt. nalpha_occupied) then 
         nszero = nbeta_occupied - nalpha_occupied 
         do a = nbeta_occupied, nalpha_occupied, -1 
            b = b + 1 
            szero(b) = a 
         enddo 
      endif 

      if (nalpha_occupied .eq. nbeta_occupied) go to 100 
c
c     Reduce second index. 
c     ------------------- 

      do i = 1,  nszero 
      do b = b1, b2  
         if (szero(i) .eq. b) then 
         do a = a1, a2 

            array2(a,b) = 0.0 

         enddo ! b 
         endif 
      enddo ! a   
      enddo ! i 
c
c     Determine range of doubly occupied indices to be removed. 
c     --------------------------------------------------------- 

      b      = 0 
      ndzero = 0 
      if (nalpha_occupied .gt. nbeta_occupied) then 
         ndzero = nbeta_occupied 
         do a = 1, nbeta_occupied  
            b = b + 1 
            dzero(b) = a 
         enddo 
      endif 
      if (nbeta_occupied .gt. nalpha_occupied) then 
         ndzero = nalpha_occupied 
         do a = 1, nalpha_occupied  
            b = b + 1 
            dzero(b) = a 
         enddo 
      endif 
c
c     Reduce first index. 
c     ------------------- 

      do i = 1,  ndzero 
      do a = a1, a2  
         if (dzero(i) .eq. a) then 
         do b = b1, b2 

            array2(a,b) = 0.0 

         enddo ! a 
         endif 
      enddo ! b   
      enddo ! i 

100   continue 
c
      return
      end

