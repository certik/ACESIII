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
      subroutine dump_block(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Print the first and last sample of a block.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
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
      integer array, index(mx_array_index), nindex, ierr
      integer block, blkndx, seg(mx_array_index)
      integer find_current_block
      integer*8 indblk, get_block_index
      integer*8 daddr, get_index_from_base
      integer type
      integer stack
      
      integer comm

      integer val1(mx_array_index), val2(mx_array_index)

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer n, nb, proc, map, company_comm
      integer na1,na2,ni1,ni2
      integer pst_get_company_comm, pst_get_company_rank
      company_comm    = pst_get_company_comm(me)
      proc  = pst_get_company_rank(me)

#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      array = op(c_result_array)
      type  = array_table(c_array_type,array)
      nb = array_table(c_numblks, array)
      map = array_table(c_block_map, array)

      nindex = array_table(c_nindex, array)
      n = 1
      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1,array)
         seg(i) = index_table(c_current_seg,index(i))

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

         call get_index_segment(index(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
         n = n * (val2(i)-val1(i)+1)
      enddo

c---------------------------------------------------------------------------
c   Get array data address.
c---------------------------------------------------------------------------

      if (type .eq. static_array) then
         daddr = address_table(array)
         indblk = get_index_from_base(daddr, x, 2)
         print *,'Dump of static array ',array,' daddr = ',daddr,
     *      ' range = ',(val1(i),val2(i),i=1,nindex)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)
         print *,'Task ',me,' at line number ',current_line,
     *     ' ARRAY ',array,' BLOCK ',block,' ind = ',
     *   (index(i),i=1,nindex),' seg = ',(seg(i),i=1,nindex),
     *   ' actual size = ',n 
         if (block .lt. 0) return 

         if (type .eq. distributed_array .and.
     *       block_map_table(c_processor,map+block-1).eq.proc) then
            call dump_block_data(array,block)
         else
            call dump_block_data(array,block)
         endif

         stack = array_table(c_array_stack, array)
         indblk = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif

      if (nindex .eq. 4) then
         call print_samples(me, array, block, 
     *                   x(indblk), val1(1),val2(1),
     *                   val1(2),val2(2),
     *                   val1(3),val2(3),
     *                   val1(4),val2(4))
      else if (nindex .eq. 2) then
         if (type .eq. static_array) then
            na1 = array_table(c_index_range1,array)
            na2 = array_table(c_index_range2,array)
            ni1 = array_table(c_index_range1+1,array)
            ni2 = array_table(c_index_range2+1,array)
         else
            na1 = val1(1)
            na2 = val2(1)
            ni1 = val1(2)
            ni2 = val2(2)
         endif
         call print_samples2(me, array, block,
     *                   x(indblk), val1(1),val2(1),
     *                   val1(2),val2(2), na1,na2,ni1,ni2)
      endif
      
      return
      end

      subroutine print_samples(me, array, block, x, 
     *                           a1,a2,i1,i2,b1,b2,j1,j2)
      implicit none
      include 'trace.h'

      integer me, array, block
      integer a1,a2,b1,b2,i1,i2,j1,j2
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)


      integer a,b,i,j
      double precision eps
      double precision val
      double precision xmax, xmin, sumsq

      xmax = x(a1,i1,b1,j1)
      xmin = x(a1,i1,b1,j1)
      sumsq = 0.
      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         xmax = max(x(a,i,b,j), xmax)
         xmin = min(x(a,i,b,j), xmin)
         sumsq = sumsq + x(a,i,b,j)*x(a,i,b,j)
      enddo
      enddo
      enddo
      enddo
 
      print *,'Task ',me,' line ',current_line,
     *        ' Array, block = ',array,block,
     *        ' first = ',
     *        x(a1,i1,b1,j1), 
     *        ' last = ',x(a2,i2,b2,j2),
     *        ' xmin, xmax, sumsq = ',xmin,xmax,sumsq,
     *        ' ranges ',a1,a2,i1,i2,b1,b2,j1,j2
      return
      end

      subroutine print_samples2(me, array, block, x, 
     *                           a1,a2,i1,i2, na1, na2,ni1,ni2)
      implicit none
      include 'trace.h'

      integer me, array, block
      integer a1,a2,i1,i2
      integer na1,na2,ni1,ni2
      double precision x(na1:na2,ni1:ni2)


      integer a,i
      double precision eps
      double precision val
      double precision xmax, xmin, sumsq

      xmax = x(a1,i1)
      xmin = x(a1,i1)
      sumsq = 0.
      do i = i1,i2
      do a = a1,a2
         xmax = max(x(a,i), xmax)
         xmin = min(x(a,i), xmin)
         sumsq = sumsq + x(a,i)*x(a,i)
      enddo
      enddo
 
      print *,'Task ',me,' line ',current_line,
     *        ' Array, block = ',array,block,
     *        ' first = ',
     *        x(a1,i1), 
     *        ' last = ',x(a2,i2),
     *        ' xmin, xmax, sumsq = ',xmin,xmax,sumsq,
     *        ' ranges ',a1,a2,i1,i2
c      do i = i1,i2
c      do a = a1,a2
c         print *,'Task ',me,' array, a,b ',array,a,i,' value ',
c     *         x(a,i)
c      enddo
c      enddo
      return
      end

