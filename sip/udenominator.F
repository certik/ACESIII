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
      subroutine udenominator(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Divides each sample of the block given by the array argument by the 
c   denominator
c       eps = epsilon(p)-epsilon(q)
c   The p=q value is set to zero.  
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

      integer i, j, k
      integer array, index, nindex, ierr
      integer block, blkndx, seg
      integer find_current_block
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

      if (dryrun) return
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

      if (nindex .eq. 4) then
         call epsu_divide(x(indblk), val1(1),val2(1),type(1),
     *                    val1(2),val2(2),type(2),
     *                    val1(3),val2(3),type(3),
     *                    val1(4),val2(4),type(4), 
     *                    epsilon,epsilonb)
      else if (nindex .eq. 2) then
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

         call epsu_divide2(x(indblk), val1(1),val2(1),type(1),
     *                     na1, na2,
     *                     val1(2),val2(2),type(2),ni1,ni2,
     *                     epsilon,epsilonb)    
      else
         print *,'Error in energy_denominator: Result array ',
     *      'must have either 4 or 2 indices.'
         print *,'array ',array,' has ',nindex,' indices.'
         call abort_job()
      endif
      
      return
      end

      subroutine epsu_divide(x, a1,a2,atype,
     *                       i1,i2,itype,b1,b2,btype,
     *                       j1,j2, jtype, epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j
      double precision eps, epsa, epsb, epsi, epsj
      double precision val

      do j = j1,j2
         if (jtype .eq. mobindex) then
            epsj = epsilonb(j)
         else
            epsj = epsilon(j)
         endif
      do b = b1,b2
         if (btype .eq. mobindex) then
            epsb = epsilonb(b)
         else
            epsb = epsilon(b)
         endif
      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)
         else
            epsi = epsilon(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif
         val =  x(a,i,b,j)
         eps = epsi + epsj - epsa - epsb
         if (eps .ne. 0.0) x(a,i,b,j) = val/eps

c A more robust check may be needed here VFL 

      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine epsu_divide2(x, a1,a2,atype,na1,na2,i1,i2,itype,
     *                        ni1,ni2,epsilon, epsilonb)
      implicit none
      include 'interpreter.h'

      integer a1,a2,i1,i2
      integer atype, itype
      integer na1, na2, ni1, ni2
      double precision x(na1:na2,ni1:ni2)

      double precision epsilon(*), epsilonb(*)

      integer a,i
      double precision eps, epsa, epsi
      double precision val

      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)
         else
            epsi = epsilon(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            epsa = epsilonb(a)
         else
            epsa = epsilon(a)
         endif
         val =  x(a,i)
         eps = epsi - epsa

c        if (eps .ne. 0.0) x(a,i) = val/eps

         if (dabs(eps) .gt. 1.0d-10) x(a,i) = val/eps  
         if (dabs(eps) .lt. 1.0d-10) then 
            x(a,i) = 0.0d0  
            if (a .ne. i) write(6,*) ' Small denominator may need  
     *         to be eliminated in energy_denominator', 
     *         'a, i', a,i, 'den' ,eps, 'value', x(a,i) 
            x(a,i) = 0.0d0  
         endif 

c If indices are the same do not modify array 

         if (a .eq. i) x(a,i) = 0.0d0  
c
      enddo
      enddo

      return
      end
