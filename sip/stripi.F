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
      subroutine stripi(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c--------------------------------------------------------------------------
c  In this routine the index i is stripped from the array1 and the 
c  resulting array is returned as array2.   
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
      integer*8 indblk1, indblk2, get_block_index
      integer stack
      integer nsend 

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index) 
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
c   Detrmine the parameters of the first array: c_result_array 
c---------------------------------------------------------------------------

      array = op(c_result_array)
      nindex = array_table(c_nindex, array) 

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         sindex(i) = index 
         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, fop1(i), fop2(i))
      enddo


      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk1 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk1 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif

c---------------------------------------------------------------------------
c   Detrmine the parameters of the second array: c_op1_array 
c---------------------------------------------------------------------------

      array = op(c_op1_array)
      nindex = array_table(c_nindex, array) 

      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)

         findex(i) = index 
         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, sop1(i), sop2(i))
      enddo


      if (array_table(c_array_type,array) .eq. static_array) then
         addr    = address_table(array)
         indblk2 = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack   = array_table(c_array_stack,array)
         indblk2 = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif

c     write(6,*) ' ******************** ' 
c     write(6,*) ' OP1 :', (sindex(i), i=1,4) 
c     write(6,*) ' OP2 :', (findex(i), i=1,4) 

200   continue 

      nsend = 1 
      do i = 1, nindex  
         nsend = nsend*(sop2(i)-sop1(i)+1) 
      enddo 
      if (nindex .eq. 4) then 
      call do_istrip4(x(indblk1), sindex, 
     *                  fop1(1), fop2(1), fop1(2), fop2(2), 
     *                  fop1(3), fop2(3), fop1(4), fop2(4),
     *                  x(indblk2), findex,  
     *                  sop1(1), sop2(1), sop1(2), sop2(2), 
     *                  sop1(3), sop2(3), sop1(4), sop2(4), 
     *                  nsend)
      endif 

      if (nindex .eq. 3) then 
      call do_istrip3(x(indblk1), sindex, 
     *                  fop1(1), fop2(1), fop1(2), fop2(2), 
     *                  fop1(3), fop2(3),
     *                  x(indblk2), findex,  
     *                  sop1(1), sop2(1), sop1(2), sop2(2), 
     *                  sop1(3), sop2(3), 
     *                  nsend)
      endif 
 
      return
      end

      subroutine do_istrip4(x, sindex, a1, a2, b1, b2, c1, c2, d1, d2,
     *                      y, findex, e1, e2, f1, f2, g1, g2, h1, h2, 
     *                      nsend)
c--------------------------------------------------------------------------
c  Routine which performs the actual stripping. Only ONE index can be 
c  stripped and the arrays must be 4-D.   
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer e1,e2,f1,f2,g1,g2,h1,h2
      integer sindex(mx_array_index), findex(mx_array_index), 
     *        match(4) 

      integer i, j, k, a, b, c, d 
      integer m, n, lda, ldb, ldc 
      integer nsend, strip  

      double precision x(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision y(e1:e2,f1:f2,g1:g2,h1:h2)

      nsend = (e2-e1+1)*(f2-f1+1)*(g2-g1+1)*(h2-h1+1) 

c *********************************************************************************** 
c     Find the nonmatching indeces and make sure that array2 has simple indeces 
c     in those locations. 
c *********************************************************************************** 

      do i = 1, 4 
         if (sindex(i) .ne. findex(i)) then 
             if (i .eq. 1) match(i) = e1 
             if (i .eq. 2) match(i) = f1 
             if (i .eq. 3) match(i) = g1 
             if (i .eq. 4) match(i) = h1 
             strip = i 
         else 
             match(i) = 0 
         endif 
      enddo 

c *********************************************************************************** 
c     The fourth index stripped 
c *********************************************************************************** 

      if ((strip .eq. 4) .and. (h1 .ge. d1) .and. (h1 .le. d2)) then 

      do c = c1, c2 
      do b = b1, b2 
      do a = a1, a2 
c     do d = d1, d2 

         y(a,b,c,h1) = x(a,b,c,h1) 

c     enddo 
      enddo 
      enddo 
      enddo 

      go to 100 

      endif 

c *********************************************************************************** 
c     The Third index stripped 
c *********************************************************************************** 

      if ((strip .eq. 3) .and. (g1 .ge. c1) .and. (g1 .le. c2)) then 

c     do c = c1, c2 
      do d = d1, d2 
      do b = b1, b2 
      do a = a1, a2 

         y(a,b,g1,d) = x(a,b,g1,d) 

      enddo 
      enddo 
      enddo 
c     enddo 

      go to 100 

      endif 

c *********************************************************************************** 
c     The Second index stripped 
c *********************************************************************************** 

      if ((strip .eq. 2) .and. (f1 .ge. b1) .and. (f1 .le. b2)) then 

      do d = d1, d2 
      do c = c1, c2 
c     do b = b1, b2 
      do a = a1, a2 

         y(a,f1,c,d) = x(a,f1,c,d) 

      enddo 
      enddo 
c     enddo 
      enddo 

      go to 100 

      endif 

c *********************************************************************************** 
c     The First index stripped 
c *********************************************************************************** 

      if ((strip .eq. 1) .and. (e1 .ge. a1) .and. (e1 .le. a2)) then 

      do d = d1, d2 
      do c = c1, c2 
      do b = b1, b2 
c     do a = a1, a2 

         y(e1,b,c,d) = x(e1,b,c,d) 

      enddo 
      enddo 
c     enddo 
      enddo 

      endif 

      write(6,*) ' No indeces were found to strip in strip4 '
      write(6,*) ' Strip = ', strip  
      write(6,*) ' Match(i) = ', (match(i), i=1,4)  
      call abort_job() 

100   continue 

      return
      end


      subroutine do_istrip3(x, sindex, a1, a2, b1, b2, c1, c2,
     *                      y, findex, e1, e2, f1, f2, g1, g2,  
     *                      nsend)
c--------------------------------------------------------------------------
c  Routine which performs the actual stripping. Only ONE index can be 
c  stripped and the arrays must be 4-D.   
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer e1,e2,f1,f2,g1,g2,h1,h2
      integer sindex(mx_array_index), findex(mx_array_index), 
     *        match(mx_array_index) 

      integer i, j, k, a, b, c, d 
      integer m, n, lda, ldb, ldc 
      integer nsend, strip  

      double precision x(a1:a2,b1:b2,c1:c2)
      double precision y(e1:e2,f1:f2,g1:g2)

      nsend = (e2-e1+1)*(f2-f1+1)*(g2-g1+1) 

c *********************************************************************************** 
c     Find the nonmatching indeces and make sure that array2 has simple indeces 
c     in those locations. 
c *********************************************************************************** 

      do i = 1, 3 
         if (sindex(i) .ne. findex(i)) then 
             if (i .eq. 1) match(i) = e1 
             if (i .eq. 2) match(i) = f1 
             if (i .eq. 3) match(i) = g1 
             strip = i 
         else 
             match(i) = 0 
         endif 
      enddo 

c *********************************************************************************** 
c     The Third index stripped 
c *********************************************************************************** 

      if ((strip .eq. 3) .and. (g1 .ge. c1) .and. (g1 .le. c2)) then 

c     do c = c1, c2 
      do b = b1, b2 
      do a = a1, a2 

         y(a,b,g1) = x(a,b,g1) 

      enddo 
      enddo 
c     enddo 

      go to 100 

      endif 

c *********************************************************************************** 
c     The Second index stripped 
c *********************************************************************************** 

      if ((strip .eq. 2) .and. (f1 .ge. b1) .and. (f1 .le. b2)) then 

      do c = c1, c2 
c     do b = b1, b2 
      do a = a1, a2 

         y(a,f1,c) = x(a,f1,c) 

      enddo 
      enddo 
c     enddo 

      go to 100 

      endif 

c *********************************************************************************** 
c     The First index stripped 
c *********************************************************************************** 

      if ((strip .eq. 1) .and. (e1 .ge. a1) .and. (e1 .le. a2)) then 

      do c = c1, c2 
      do b = b1, b2 
c     do a = a1, a2 

         y(e1,b,c) = x(e1,b,c) 

      enddo 
      enddo 
c     enddo 

      go to 100 

      endif 

      write(6,*) ' No indeces were found to strip in strip3 '
      write(6,*) ' Strip = ', strip  
      write(6,*) ' Match(i) = ', (match(i), i=1,3)  
      call abort_job() 

100   continue 

      return
      end

