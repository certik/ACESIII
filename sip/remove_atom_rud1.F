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
      subroutine remove_atom_rud1(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c--------------------------------------------------------------------------
c Takes a four-dimensional array array(a,b,c,d) and removes the elements
c in which atom(b) = atom(a,c,or d).    
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
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

         call doremove_atom_rud1(x(indblk), val1(1),val2(1),type(1),
     *                   val1(2),val2(2),type(2),
     *                   val1(3),val2(3),type(3),
     *                   val1(4),val2(4),type(4))
      else
         print *,'Error in remove_atom_rud1: Result array ',
     *      'must have either 4 indices.'
         print *,'array ',array,' has ',nindex,' indices.'
         call abort_job()
      endif

      return
      end

      subroutine doremove_atom_rud1(x, a1,a2,atype,
     *                      i1,i2,itype,b1,b2,btype,
     *                      j1,j2, jtype)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h' 

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      integer a,b,c,d,i,j, m1, m2, n1, n2, r1, r2, s1, s2 
      integer m, n, r, s, aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2 
      double precision val

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through to determine
c   the atomic labels .
c-----------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, i1, n1)
      call lookup_shell(end_nfps, nshells, i2, n2)
      call lookup_shell(end_nfps, nshells, b1, r1)
      call lookup_shell(end_nfps, nshells, b2, r2)
      call lookup_shell(end_nfps, nshells, j1, s1)
      call lookup_shell(end_nfps, nshells, j2, s2)

      do m = m1, m2 
            if (m .eq. 1) then
               aa1 = 1
            else
               aa1 = end_nfps(m-1) + 1
            endif
            aa2 = end_nfps(m) 
         do n = n1, n2
            if (n .eq. 1) then
               bb1 = 1
            else
               bb1 = end_nfps(n-1) + 1
            endif
            bb2 = end_nfps(n)
         do r = r1, r2
            if (r .eq. 1) then
               cc1 = 1
            else
               cc1 = end_nfps(r-1) + 1
            endif
            cc2 = end_nfps(r)
         do s = s1, s2
               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)  

             if ((atom(n) .eq. atom(m)) .or. (atom(n) .eq. atom(r)) .or. 
     *           (atom(n) .eq. atom(s))) then 

             do a = aa1, aa2  
             do b = bb1, bb2  
             do c = cc1, cc2  
             do d = dd1, dd2  
                x(a,b,c,d) = 0.0d0 
             enddo 
             enddo 
             enddo 
             enddo 

             endif 

         enddo ! m 
         enddo ! n 
         enddo ! r 
      enddo ! s 


      return
      end


