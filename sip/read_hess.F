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
      subroutine read_hess(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c     The 4-index array (hessian at present) is read in and summed 
c     into an array in a common block. 
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
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

      integer i, j, k
      integer array, index, nindex, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer*8 indblk, get_block_index
      integer stack
      
      integer comm

      integer val1(mx_array_index), val2(mx_array_index)
      integer*8 addr, get_index_from_base

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
      write_hess = .true.

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      array = op(c_result_array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         index = array_table(c_index_array1+i-1,array)
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

      if (nindex .ne. 4) then
         print *,'Error in read_hess: ',
     *      'Hessian must be a 4-index array!.'
         print *,'array ',array,' has ',nindex,' indices.'
         call abort_job()
      endif

      call add_hess(x(indblk),val1(1),val2(1),val1(2),val2(2),
     *                        val1(3),val2(3),val1(4),val2(4),
     *                        x(ix_hess))
      
      return
      end


      subroutine add_hess(x,iatoms,iatome,icomponents,icomponente,
     *                      jatoms,jatome,jcomponents,jcomponente, 
     *                      hess) 
      implicit none
      include 'parallel_info.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'

      integer i,j,ihess,jhess 
      integer iatoms,iatome,icomponents,icomponente  
      integer jatoms,jatome,jcomponents,jcomponente  
      double precision x(iatoms:iatome,icomponents:icomponente,
     *                   jatoms:jatome,jcomponents:jcomponente)
      double precision hess(3*ncenters,3*ncenters) 


      if (dbg) then
         write(6,*) ' HCOORD :', iatoms,icomponents,jatoms,
     *                    jcomponents
         write(6,*) ' Hessian from Loops ',x(iatoms,icomponents, 
     *                                    jatoms,jcomponents)  
      endif

      if (me .ne. 0) return  

      ihess = (iatoms-1)*3 + icomponents 
      jhess = (jatoms-1)*3 + jcomponents 

      hess(jhess,ihess) = hess(jhess,ihess) 
     *                  + x(iatoms,icomponents,jatoms,jcomponents) 

c     if (ihess .ne. jhess) 
c      hess(jhess,ihess) = hess(jhess,ihess) 
c     *                  + x(iatoms,icomponents,jatoms,jcomponents) 

33    format(3F16.8)

      return
      end


