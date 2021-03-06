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
      subroutine scont1(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Performs the contraction of derivative 1-electron integrals in the AO
c   domain with an intermediate term derived from the back-transformed
c   gamma matrix. The contraction is performed on an atom-by-atom basis
c   and summed into the appropriate component of the gradient array.
c
c   scont1 must be executed with the following syntax:
c   execute scont1 i(mu,nu)
c
c   The block i(mu,nu) is the intermediate term which is 
c   to be contracted with the "S" derivative integrals for its 
c   block of (mu,nu) across all the atoms in the molecule.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'dbugcom.h'
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

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
      integer array, nind, ind(mx_array_index), seg(mx_array_index)
      integer a1, a2, b1, b2
      integer i, j, m, n
      integer shandle, iatom
      integer*8 iscr, indblk
      integer find_current_block
      integer*8 get_block_data_index
      integer block, blkndx, stack
      integer igrad, component
      integer allocate_scratch_block
      integer ierr

      double precision dummy, y
      integer flopcount
      integer a,b

#ifdef ALTIX
      dptr = dshptr
#endif

      if (compute_1e_integrals) then
         compute_1e_integrals = .false.

         if (me .eq. 0) then

c-------------------------------------------------------------------------
c   Sum the nuclear-nuclear repulsion term into the gradient.
c--------------------------------------------------------------------------

            iatom = 0  
            j     = 0
            do i = 1, 3*ncenters
               j = j + 1 
               if (j .gt. 3) j = 1 
               if (j .eq. 1) iatom = iatom + 1
               gradient_data(i) = gradient_data(i) + 
     *                            NNgrad(j,iatom)
            enddo

            if (dbg) then
               print *,'SCONT1 on Task 0: Summed NNgrad into Gradient'
               do i = 1, 3*ncenters
                  print *,'Gradient i: ',i,gradient_data(i)
               enddo
            endif
         endif
      endif

c---------------------------------------------------------------------------
c   Allocate a scratch block to use in the contractions.
c---------------------------------------------------------------------------

      ierr =  allocate_scratch_block(x, iscr, shandle, array_table,
     *                            narray_table, index_table, 
     *                            nindex_table, block_map_table)
      if (ierr .ne. 0) then
         print *,'Error in SCONT1: Cannot allocate scratch block'
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Find the ranges of the AO indices of the input block.
c--------------------------------------------------------------------------

         array = op(c_result_array)

         nind  = array_table(c_nindex,array)
         if (nind .ne. 2) then
            print *,'Error: SCONT1 requires a 2-index argument array.'
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
         if (block .le. 0) then
            print *,'Task ',me,
     *         ' Error: SCONT1 cannot find current block',
     *         ' of array ',array,' indices ',ind(1),ind(2),' segs ',
     *         seg(1),seg(2)
            call abort_job()
         endif

         stack = array_table(c_array_stack, array)
         indblk = get_block_data_index(array, block, stack,
     *                                        blkndx, x)

c---------------------------------------------------------------------------
c   Perform contraction of the data block with the integrals.
c---------------------------------------------------------------------------

      call contract_oed_derv_integrals(3, x(indblk), x(iscr),
     *                  a1, a2, b1, b2)

      call free_scratch_block(shandle)
      return
      end

