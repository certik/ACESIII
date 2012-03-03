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
      subroutine compute_Sderivative_integrals(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   CHANGE COMMENT 
c   Multiplies the each sample of the block given by the array argument by 
c   diagonal eigenvalue matrix.  The eigenvalues matrix is not stored in
c   a 2-dimensional form, since it can be represented by a vector of
c   the diagonal elements.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'hess.h'
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

      integer array, ind(mx_array_index), nindex, ierr
      integer block, blkndx, shandle, seg
      integer*8 iderblk(4,3), indblk, get_block_data_index
      integer stack, handle, handle2, icenter, match
      integer find_current_block

      integer val1(mx_array_index), val2(mx_array_index)
      integer iatom, igrad, flopcount
      integer component

      integer i, j, k
      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer index 
      integer nsend
      integer der_flags(12)
      
      integer my_comm_rank, comm
      integer*8 c_loc64, ixx

      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (iptr, ibuf)
#else
      equivalence (ibuf, dbuf)
      common dbuf
#endif

      integer flag, iarray 

      if (dryrun) return
#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c----------------------------------------------------------------------------
c   Locate the data for the intermediate block.
c---------------------------------------------------------------------------

      array = op(c_result_array)
      nindex = array_table(c_nindex, array)
      if (nindex .ne. 4) then 
         write(6,*) ' ERROR : number of indices wrong when calling ',
     *              ' compute_Sderivative_integrals ', nindex 
         call abort_job()
      endif 

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
         seg = index_table(c_current_seg,ind(i))
         call get_index_segment(ind(i), seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

      a1 = val1(1)
      a2 = val2(1)
      b1 = val1(2)
      b2 = val2(2)
      c1 = val1(3)
      c2 = val2(3)
      d1 = val1(4)
      d2 = val2(4)

      block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

      stack = array_table(c_array_stack,array)
      indblk = get_block_data_index(array, block, stack,
     *                                        blkndx, dbuf)

c-----------------------------------------------------------------------------
c   Calculate the integral blocks.
c-----------------------------------------------------------------------------

         do k = 1, 12
            der_flags(k) = 0
         enddo

c-----------------------------------------------------------------------------
c   Computation of Hessian contribution.
c-----------------------------------------------------------------------------

         call contract_erd_2der(a1,a2,b1,b2,c1,c2,d1,d2,
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               dbuf(indblk), dbuf(ix_hess))

      go to 999
c-----------------------------------------------------------------------------
c   Test of second-derivatives.
c-----------------------------------------------------------------------------

         write(6,*) ' Computing Second-derivatives '
         call tester_MRT2(a1,a2,b1,b2,c1,c2,d1,d2,
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               dbuf(indblk), dbuf(ix_hess))

c-----------------------------------------------------------------------------
c   Test of third-derivatives.
c-----------------------------------------------------------------------------

         write(6,*) ' Computing Third-derivatives '
         call tester_MRT3(a1,a2,b1,b2,c1,c2,d1,d2,
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               dbuf(indblk))
c-----------------------------------------------------------------------------
c   Test of fourth-derivatives.
c-----------------------------------------------------------------------------

         write(6,*) ' Computing Fourth-derivatives '
         call tester_MRT4(a1,a2,b1,b2,c1,c2,d1,d2,
     *               der_flags,
     *               dbuf(dpkgscr),zmax, ibuf(ipkgscr),
     *               dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *               ibuf(iccbeg), ibuf(iccend),
     *               dbuf(indblk))

999   continue

      return
      end
