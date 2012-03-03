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
      subroutine hcontxy(array_table, 
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
c   hcontxy must be executed with the following syntax:
c   execute hcontxy i(mu,nu)
c
c   The block i(mu,nu) is the intermediate term which is 
c   to be contracted with the second derivative "H" integrals for its 
c   block of (mu,nu) across all the atoms in the molecule.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'hess.h'
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
      integer i, m, n, nl2
      integer shandle, iatom
      integer*8 iscr, indblk
      integer find_current_block
      integer*8 get_block_data_index
      integer block, blkndx, stack
      integer igrad, component, intermediate
      integer allocate_scratch_block, ierr

      double precision dummy, y
      integer flopcount

      nl2 = (nbasis * nbasis + nbasis)/2
      write_hess = .true.

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Allocate a scratch block to use in the contractions.
c---------------------------------------------------------------------------

      ierr = allocate_scratch_block(x, iscr, shandle, array_table,
     *                            narray_table, index_table, 
     *                            nindex_table, block_map_table)
      if (ierr .ne. 0) then
         print *,'Error in HCONT1; Cannot allocate scratch block'
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Find the ranges of the AO indices of the input block.
c--------------------------------------------------------------------------

         array = op(c_result_array)

         nind  = array_table(c_nindex,array)
         if (nind .ne. 2) then
            print *,'Error: HCONT1 requires a 2-index argument array.'
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
c   Perform contraction of the data block with the integrals.
c---------------------------------------------------------------------------

      if (intpkg .eq.flocke_package) then


c        call contract_oed_2derv_type1a(1, x(indblk), x(iscr),
c    *                  a1, a2, b1, b2, x(ix_hess))

         call contract_oed_2derv_type1b(1, x(indblk), x(iscr),
     *                  a1, a2, b1, b2, x(ix_hess))
        call contract_oed_2derv_integrals(2, x(indblk), x(iscr),
     *                  a1, a2, b1, b2, x(ix_hess))

c         call contract_oed_2derv_numint(2, x(indblk), x(iscr),
c     *                  a1, a2, b1, b2, x(ix_hess))

      endif

      call free_scratch_block(shandle)
      return
      end
