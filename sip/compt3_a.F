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
      subroutine compt3_a(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c  The (T 4v) calculation is performed.  
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'machine_types.h'
      include 'saved_data.h'
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
      integer array, array_type, nscalar_table, sind 
      integer index, nindex, seg 
      double precision scalar_table(nscalar_table), tind 
      integer*8 address_table(narray_table)
      integer*8 addr, get_index_from_base 
      integer*8 indblk, indblka(20), get_block_data_index, 
     *          get_block_index
      integer type(mx_array_index)

      integer Top1(mx_array_index), Top2(mx_array_index), 
     *        Tindex(mx_array_index),   
     *        Sindex(18,mx_array_index)  
      integer Sop1(20,mx_array_index), Sop2(20,mx_array_index) 
      integer S1, S2 

      integer Na, nsend, i, j, k
      integer flag, handle, ierr
      integer blk, blkndx, stack, find_current_block, block  

      integer shandle(20), istat
      integer allocate_scratch_block
      integer*8 iscr(20)
      integer*8 c_loc64 !, n

      double precision x(1)
      double precision y(1)
      double precision scr(1)
      integer t3ablocks  
      common /T3a_blocks/Na, t3ablocks(20) 


c     save shandle 
c     save iscr 

#ifdef ALTIX
      pointer (dptr, x)
      pointer (dptr, y)
#else
      common x
      common y
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

c------------------------------------------------------------------------------
c  Determine the parameters of the output array T3(a,i,b,j,c,k): c_result_array
c------------------------------------------------------------------------------
      
      array  = op(c_result_array)
      nindex = array_table(c_nindex, array)

      if (nindex .ne. 6) then 
         write(6,*) ' Output from Ta must have six indeces, you have:', 
     *                nindex 
         call abort_job() 
      endif 
      
      do i = 1, nindex
         index   = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg     = index_table(c_current_seg,index)
      
         Tindex(i) = index
         call get_index_segment(index, seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, Top1(i), Top2(i))
      enddo
      
      block = find_current_block(array, array_table(1,array),
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           block_map_table, blkndx)
      
      stack  = array_table(c_array_stack,array)
      indblk = get_block_index(array, block, stack,
     *                         blkndx, y, .true.)

c------------------------------------------------------------------------------
c  Determine the parameters of the input/contracted arrays.
c------------------------------------------------------------------------------

      do i = 1, 12 
         array = t3ablocks(i) 
         nindex = array_table(c_nindex, array)

         if (nindex .ne. 4) then 
            write(6,*) ' Input arrays to Ta must have four indeces, 
     *                   you have:', nindex 
            call abort_job() 
         endif 
      
         do j = 1, nindex
            index   = array_table(c_index_array1+j-1,array)
            type(j) = index_table(c_index_type, index)
            seg     = index_table(c_current_seg,index)
      
            Sindex(i,j) = index
            call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, S1, S2)
            Sop1(i,j) = S1 
            Sop2(i,j) = S2 
         enddo

         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack = array_table(c_array_stack,array)
         indblka(i) = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      enddo

c------------------------------------------------------------------------
c   Allocate scratch blocks used to store the sorted 4-index arrays.  
c------------------------------------------------------------------------
         
       do i = 1, 15 
         istat = allocate_scratch_block(scr, iscr(i), shandle(i),
     *                           array_table, narray_table,
     *                           index_table, nindex_table,
     *                           block_map_table)
         if (istat .ne. 0) then
            print *,'Error: Unable to allocate scratch blocks for T.'
            call abort_job()
         endif
      enddo 

c------------------------------------------------------------------------
c   Perform the actual calculation.  
c------------------------------------------------------------------------

      call do_compt3_a1(y(indblk), Top1, Top2, 
     *                 x(indblka(1)),     
     *             Sop1(1,1), Sop1(1,2), Sop1(1,3), Sop1(1,4),
     *             Sop2(1,1), Sop2(1,2), Sop2(1,3), Sop2(1,4),
c
     *                 x(indblka(2)),   
     *             Sop1(2,1), Sop1(2,2), Sop1(2,3), Sop1(2,4),
     *             Sop2(2,1), Sop2(2,2), Sop2(2,3), Sop2(2,4),
c
     *                 x(indblka(3)),   
     *             Sop1(3,1), Sop1(3,2), Sop1(3,3), Sop1(3,4),
     *             Sop2(3,1), Sop2(3,2), Sop2(3,3), Sop2(3,4),
c
     *                 x(indblka(4)),   
     *             Sop1(4,1), Sop1(4,2), Sop1(4,3), Sop1(4,4),
     *             Sop2(4,1), Sop2(4,2), Sop2(4,3), Sop2(4,4),
c
     *                 x(indblka(5)),   
     *             Sop1(5,1), Sop1(5,2), Sop1(5,3), Sop1(5,4),
     *             Sop2(5,1), Sop2(5,2), Sop2(5,3), Sop2(5,4),
c
     *                 x(indblka(6)),   
     *             Sop1(6,1), Sop1(6,2), Sop1(6,3), Sop1(6,4),
     *             Sop2(6,1), Sop2(6,2), Sop2(6,3), Sop2(6,4),
c
     *                 x(indblka(7)),   
     *             Sop1(7,1), Sop1(7,2), Sop1(7,3), Sop1(7,4),
     *             Sop2(7,1), Sop2(7,2), Sop2(7,3), Sop2(7,4),
c
     *                 x(indblka(8)),   
     *             Sop1(8,1), Sop1(8,2), Sop1(8,3), Sop1(8,4),
     *             Sop2(8,1), Sop2(8,2), Sop2(8,3), Sop2(8,4),
c
     *                 x(indblka(9)),   
     *             Sop1(9,1), Sop1(9,2), Sop1(9,3), Sop1(9,4),
     *             Sop2(9,1), Sop2(9,2), Sop2(9,3), Sop2(9,4),
c
     *                 x(indblka(10)),   
     *             Sop1(10,1), Sop1(10,2), Sop1(10,3), Sop1(10,4),
     *             Sop2(10,1), Sop2(10,2), Sop2(10,3), Sop2(10,4),
c
     *                 x(indblka(11)),   
     *             Sop1(11,1), Sop1(11,2), Sop1(11,3), Sop1(11,4),
     *             Sop2(11,1), Sop2(11,2), Sop2(11,3), Sop2(11,4),
c
     *                 x(indblka(12)),   
     *             Sop1(12,1), Sop1(12,2), Sop1(12,3), Sop1(12,4),
     *             Sop2(12,1), Sop2(12,2), Sop2(12,3), Sop2(12,4),
c
     *             Sindex, scr(iscr(1)), scr(iscr(2)), scr(iscr(3)), 
     *             scr(iscr(4)), scr(iscr(5)), scr(iscr(6)), 
     *             scr(iscr(7)), scr(iscr(8)), scr(iscr(9)), 
     *             scr(iscr(10)), scr(iscr(11)), scr(iscr(12)), 
c 
     *             nsend)

      do i = 1, 15 
         call free_scratch_block(shandle(i)) 
      enddo 

      return
      end

      subroutine do_compt3_a1(x, Top1, Top2, 
     *         LT001aa   , e1, f1, g1, h1, e1f, f1f, g1f, h1f, 
     *         LVS122aaai, e2, f2, g2, h2, e2f, f2f, g2f, h2f, 
c
     *         LT201aa   , e3, f3, g3, h3, e3f, f3f, g3f, h3f, 
     *         LVS102aaai, e4, f4, g4, h4, e4f, f4f, g4f, h4f, 
c
     *         LT101aa   , e5, f5, g5, h5, e5f, f5f, g5f, h5f, 
     *         LVS022aaai, e6, f6, g6, h6, e6f, f6f, g6f, h6f, 
c
     *         LT021aa   , e7, f7, g7, h7, e7f, f7f, g7f, h7f, 
     *         LVS120aaai, e8, f8, g8, h8, e8f, f8f, g8f, h8f, 
c
     *         LT002aa   , e9, f9, g9, h9, e9f, f9f, g9f, h9f, 
     *         LVS121aaai, e10, f10, g10, h10, e10f, f10f, g10f, h10f, 
c
     *         LT121aa   , e11, f11, g11, h11, e11f, f11f, g11f, h11f, 
     *         LVS020aaai, e12, f12, g12, h12, e12f, f12f, g12f, h12f, 
c
     *         Sindex, scr1, scr2, scr3, scr4, scr5, scr6, scr7, 
     *         scr8, scr9, scr10, scr11, scr12,  
     *         nsend)
c 
      implicit none
c     include 'interpreter.h'
c     include 'trace.h'
c     include 'mpif.h'
      include 'maxdim.h'
c     include 'saved_data.h'
c 
      integer e1, f1, g1, h1, e1f, f1f, g1f, h1f 
      integer e2, f2, g2, h2, e2f, f2f, g2f, h2f 
c
      integer e3, f3, g3, h3, e3f, f3f, g3f, h3f 
      integer e4, f4, g4, h4, e4f, f4f, g4f, h4f 
c
      integer e5, f5, g5, h5, e5f, f5f, g5f, h5f 
      integer e6, f6, g6, h6, e6f, f6f, g6f, h6f 
c
      integer e7, f7, g7, h7, e7f, f7f, g7f, h7f 
      integer e8, f8, g8, h8, e8f, f8f, g8f, h8f 
c
      integer e9, f9, g9, h9, e9f, f9f, g9f, h9f 
      integer e10, f10, g10, h10, e10f, f10f, g10f, h10f 
c
      integer e11, f11, g11, h11, e11f, f11f, g11f, h11f 
      integer e12, f12, g12, h12, e12f, f12f, g12f, h12f 
c
      integer a, b, c, i, j, k, a3, a3_s, a3_f  
      integer a1, b1, c1, i1, j1, k1  
      integer a2, b2, c2, i2, j2, k2  
      integer aa, aa1, aa2, ii, ii1, ii2 
      integer aa_range, ii_range, aa1_range, ii1_range, ii2_range, 
     *        aa2_range, a3_range   
      integer Top1(mx_array_index), Top2(mx_array_index) 
      integer ycount, ytemp, nsend,
     *        Sindex(18,mx_array_index)  
      integer y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12 
      integer y0, z0, z1, z2, z3, z4, z5, z6 
c
      double precision x(*), etemp, etemp1                      
  
      double precision LT001aa   ( e1:e1f, f1:f1f, g1:g1f, h1:h1f) 
      double precision LVS122aaai( e2:e2f, f2:f2f, g2:g2f, h2:h2f) 
c
      double precision LT201aa   ( e3:e3f, f3:f3f, g3:g3f, h3:h3f) 
      double precision LVS102aaai( e4:e4f, f4:f4f, g4:g4f, h4:h4f) 
c
      double precision LT101aa   ( e5:e5f, f5:f5f, g5:g5f, h5:h5f) 
      double precision LVS022aaai( e6:e6f, f6:f6f, g6:g6f, h6:h6f) 
c
      double precision LT021aa   ( e7:e7f, f7:f7f, g7:g7f, h7:h7f) 
      double precision LVS120aaai( e8:e8f, f8:f8f, g8:g8f, h8:h8f) 
c
      double precision LT002aa   ( e9:e9f, f9:f9f, g9:g9f, h9:h9f) 
      double precision LVS121aaai( e10:e10f,f10:f10f,g10:g10f,h10:h10f) 
c
      double precision LT121aa   ( e11:e11f,f11:f11f,g11:g11f,h11:h11f) 
      double precision LVS020aaai( e12:e12f,f12:f12f,g12:g12f,h12:h12f) 
c
      double precision scr1(*)  
      double precision scr2(*)  
      double precision scr3(*)  
      double precision scr4(*)  
      double precision scr5(*)  
      double precision scr6(*)  
      double precision scr7(*)  
      double precision scr8(*)  
      double precision scr9(*)  
      double precision scr10(*)  
      double precision scr11(*)  
      double precision scr12(*)  
c
c ------------------------------------------------------------------------- 
c     Set output array's ranges 
c ------------------------------------------------------------------------- 
c 
      a1 = Top1(1) 
      a2 = Top2(1) 
c 
      i1 = Top1(2) 
      i2 = Top2(2) 
c 
      b1 = Top1(3) 
      b2 = Top2(3) 
c 
      j1 = Top1(4) 
      j2 = Top2(4) 
c 
      c1 = Top1(5) 
      c2 = Top2(5) 
c 
      k1 = Top1(6) 
      k2 = Top2(6) 
c
c ------------------------------------------------------------------------- 
c     Determine the contracted index ranges 
c ------------------------------------------------------------------------- 
c
      do i = 1, 4 
      do j = 1, 4 
         if (Sindex(1,i) .eq. Sindex(2,j)) then 
            if (i .eq. 1) then 
               a3_s = e1  
               a3_f = e1f  
               go to 100 
            endif 
            if (i .eq. 2) then 
               a3_s = f1  
               a3_f = f1f  
               go to 100 
            endif 
            if (i .eq. 3) then 
               a3_s = g1  
               a3_f = g1f  
               go to 100 
            endif 
            if (i .eq. 4) then 
               a3_s = h1  
               a3_f = h1f  
               go to 100 
            endif 
         endif 
      enddo 
      enddo 
100   continue 
c
c ------------------------------------------------------------------------- 
c Put arrays into scratch arrays to be used sequentially 
c ------------------------------------------------------------------------- 
c
      aa_range  = a2-a1+1 
      ii_range  = i2-i1+1 
      aa1_range = b2-b1+1 
      ii1_range = j2-j1+1 
      aa2_range = c2-c1+1 
      ii2_range = k2-k1+1 
      a3_range  = a3_f-a3_s+1 

c     go to 33 
c
      ycount = 0 
      do ii1 = j1, j2 
      do ii  = i1, i2 
      do aa  = a1, a2 

         do a3 = a3_s, a3_f 
            ycount = ycount + 1 

         y1 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *      + (ii-i1)*aa_range*a3_range  
     *      + (ii1-j1)*ii_range*aa_range*a3_range 

            scr1(y1) = LT001aa(aa,ii,a3,ii1)

            if (y1 .ne. ycount) then 
            write(6,*) 'y1', y1, ycount, a3, aa, ii, ii1  
            endif 

         enddo 
      enddo 
      enddo 
      enddo 

      ycount = 0 
      do ii2 = k1, k2 
      do aa2 = c1, c2 
      do aa1 = b1, b2 
         do a3 = a3_s, a3_f 
            ycount = ycount + 1 

            y2 = (a3-a3_s+1) + (aa1-b1)*a3_range   
     *         + (aa2-c1)*aa1_range*a3_range  
     *         + (ii2-k1)*aa2_range*aa1_range*a3_range 

            if (y2 .ne. ycount) then 
            write(6,*) 'Y2', y2, ycount, a3, aa, ii, ii1  
            endif 

            scr2(y2) = LVS122aaai(aa1,a3,aa2,ii2)
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa1 = b1, b2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y4 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (aa1-b1)*aa_range*a3_range  
     *         + (ii2-k1)*aa1_range*aa_range*a3_range 
            scr4(y4) = LVS102aaai(aa1,a3,aa,ii2)

            ycount = ycount + 1 
            if (y4 .ne. ycount) then 
            write(6,*) 'Y4', y4, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do aa2 = c1, c2 
      do ii1 = j1, j2 
      do ii  = i1, i2 
         do a3 = a3_s, a3_f 
            y3 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (ii1-j1)*ii_range*a3_range  
     *         + (aa2-c1)*ii1_range*ii_range*a3_range 
            scr3(y3) = LT201aa(aa2,ii,a3,ii1)

            ycount = ycount + 1 
            if (y3 .ne. ycount) then 
            write(6,*) 'Y3', y3, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa2 = c1, c2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y6 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (aa2-c1)*aa_range*a3_range  
     *         + (ii2-k1)*aa2_range*aa_range*a3_range 

            scr6(y6) = LVS022aaai(aa,a3,aa2,ii2)

            ycount = ycount + 1 
            if (y6 .ne. ycount) then 
            write(6,*) 'Y6', y6, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii1 = j1, j2 
      do aa1 = b1, b2 
      do ii  = i1, i2 
         do a3 = a3_s, a3_f 
            y5 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa1-b1)*ii_range*a3_range  
     *         + (ii1-j1)*aa1_range*ii_range*a3_range 

            scr5(y5) = LT101aa(aa1,ii,a3,ii1)

            ycount = ycount + 1 
            if (y5 .ne. ycount) then 
            write(6,*) 'Y5', y5, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do ii1 = j1, j2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y7 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii1-j1)*aa_range*a3_range  
     *         + (ii2-k1)*ii1_range*aa_range*a3_range 

            scr7(y7) = LT021aa(aa,ii2,a3,ii1)

            ycount = ycount + 1 
            if (y7 .ne. ycount) then 
            write(6,*) 'Y7', y7, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do aa2 = c1, c2 
      do aa1 = b1, b2 
      do ii  = i1, i2 
         do a3 = a3_s, a3_f 
            y8 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa1-b1)*ii_range*a3_range  
     *         + (aa2-c1)*aa1_range*ii_range*a3_range 

            scr8(y8) = LVS120aaai(aa1,a3,aa2,ii)

            ycount = ycount + 1 
            if (y8 .ne. ycount) then 
            write(6,*) 'Y8', y8, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do ii  = i1, i2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y9 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii-i1)*aa_range*a3_range  
     *         + (ii2-k1)*ii_range*aa_range*a3_range 

            scr9(y9) = LT002aa(aa,ii,a3,ii2)

            ycount = ycount + 1 
            if (y9 .ne. ycount) then 
            write(6,*) 'Y9', y9, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do aa2 = c1, c2 
      do ii1 = j1, j2 
      do aa1 = b1, b2 
         do a3 = a3_s, a3_f 
            y10 = (a3-a3_s+1) + (aa1-b1)*a3_range   
     *         + (ii1-j1)*aa1_range*a3_range  
     *         + (aa2-c1)*ii1_range*aa1_range*a3_range 

            scr10(y10) = LVS121aaai(aa1,a3,aa2,ii1)

            ycount = ycount + 1 
            if (y10 .ne. ycount) then 
            write(6,*) 'Y10', y10, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do aa2 = c1, c2 
      do ii  = i1, i2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y12 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii-i1)*aa_range*a3_range  
     *         + (aa2-c1)*ii_range*aa_range*a3_range 
c
             scr12(y12) = LVS020aaai(aa,a3,aa2,ii)

            ycount = ycount + 1 
            if (y12 .ne. ycount) then 
            write(6,*) 'Y12', y12, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do ii1 = j1, j2 
      do aa1 = b1, b2 
         do a3 = a3_s, a3_f 
            y11 = (a3-a3_s+1) + (aa1-b1)*a3_range   
     *         + (ii1-j1)*aa1_range*a3_range  
     *         + (ii2-k1)*ii1_range*aa1_range*a3_range 
c
            scr11(y11) = LT121aa(aa1,ii2,a3,ii1)

            ycount = ycount + 1 
            if (y11 .ne. ycount) then 
            write(6,*) 'Y11', y11, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 

33    continue 
c
c ------------------------------------------------------------------------- 
c Perform calculation  
c ------------------------------------------------------------------------- 
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa2 = c1, c2 
      do ii1 = j1, j2 
      do aa1 = b1, b2 
      do ii  = i1, i2 
      do aa  = a1, a2 

         etemp  = 0 
         ycount = ycount + 1 
c
         do a3 = a3_s, a3_f 
c
            y1 = (a3-a3_s+1) 
     *         + (aa-a1)*a3_range   
     *         + (ii-i1)*aa_range*a3_range  
     *         + (ii1-j1)*ii_range*aa_range*a3_range 

            y2 = (a3-a3_s+1) + (aa1-b1)*a3_range   
     *         + (aa2-c1)*aa1_range*a3_range  
     *         + (ii2-k1)*aa2_range*aa1_range*a3_range 

            y4 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (aa1-b1)*aa_range*a3_range  
     *         + (ii2-k1)*aa1_range*aa_range*a3_range 

            y3 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (ii1-j1)*ii_range*a3_range  
     *         + (aa2-c1)*ii1_range*ii_range*a3_range 

            y6 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (aa2-c1)*aa_range*a3_range  
     *         + (ii2-k1)*aa2_range*aa_range*a3_range 

            y5 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa1-b1)*ii_range*a3_range  
     *         + (ii1-j1)*aa1_range*ii_range*a3_range 

            y7 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii1-j1)*aa_range*a3_range  
     *         + (ii2-k1)*ii1_range*aa_range*a3_range 

            y8 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa1-b1)*ii_range*a3_range  
     *         + (aa2-c1)*aa1_range*ii_range*a3_range 

            y9 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii-i1)*aa_range*a3_range  
     *         + (ii2-k1)*ii_range*aa_range*a3_range 

            y10 = (a3-a3_s+1) + (aa1-b1)*a3_range   
     *         + (ii1-j1)*aa1_range*a3_range  
     *         + (aa2-c1)*ii1_range*aa1_range*a3_range 

            y12 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii-i1)*aa_range*a3_range  
     *         + (aa2-c1)*ii_range*aa_range*a3_range 
c
            y11 = (a3-a3_s+1) + (aa1-b1)*a3_range   
     *         + (ii1-j1)*aa1_range*a3_range  
     *         + (ii2-k1)*ii1_range*aa1_range*a3_range 

            etemp = etemp 
c    *        + LT001aa(aa,ii,a3,ii1)*LVS122aaai(aa1,a3,aa2,ii2)
     *        + scr1(y1)*scr2(y2)

c    *        - LVS102aaai(aa1,a3,aa,ii2)*LT201aa(aa2,ii,a3,ii1)
     *        - scr3(y3)*scr4(y4)

c    *        - LVS022aaai(aa,a3,aa2,ii2)*LT101aa(aa1,ii,a3,ii1)
     *        - scr5(y5)*scr6(y6)

c    *        - LT021aa(aa,ii2,a3,ii1)*LVS120aaai(aa1,a3,aa2,ii)
     *        - scr7(y7)*scr8(y8)

c    *        - LT002aa(aa,ii,a3,ii2)*LVS121aaai(aa1,a3,aa2,ii1)
     *        - scr9(y9)*scr10(y10)

c    *        + LVS020aaai(aa,a3,aa2,ii)*LT121aa(aa1,ii2,a3,ii1)
     *        + scr11(y11)*scr12(y12)
c
         enddo 

         x(ycount) = etemp   

      enddo 
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 

      nsend = aa_range*ii_range*aa1_range*ii1_range*ii2_range*aa2_range    
c
c       etemp = etemp + LT001aa(aa,ii,a3,ii1)*LVS122aaai(aa1,a3,aa2,ii2)
c       etemp = etemp - LVS102aaai(aa1,a3,aa,ii2)*LT201aa(aa2,ii,a3,ii1)
c       etemp = etemp - LVS022aaai(aa,a3,aa2,ii2)*LT101aa(aa1,ii,a3,ii1)
c       etemp = etemp - LT021aa(aa,ii2,a3,ii1)*LVS120aaai(aa1,a3,aa2,ii)
c       etemp = etemp - LT002aa(aa,ii,a3,ii2)*LVS121aaai(aa1,a3,aa2,ii1)
c       etemp = etemp + LVS020aaai(aa,a3,aa2,ii)*LT121aa(aa1,ii2,a3,ii1)
c       etemp = etemp + LVS021aaai(aa,a3,aa2,ii1)*LT102aa(aa1,ii,a3,ii2)
c       etemp = etemp + LVS100aaai(aa1,a3,aa,ii)*LT221aa(aa2,ii2,a3,ii1)
c       etemp = etemp + LVS101aaai(aa1,a3,aa,ii1)*LT202aa(aa2,ii,a3,ii2)

      return
      end

      subroutine do_compt3_a3(x, Top1, Top2, 
     *         LT102aa   , e13, f13, g13, h13, e13f, f13f, g13f, h13f, 
     *         LVS021aaai, e14, f14, g14, h14, e14f, f14f, g14f, h14f, 
c
     *         LT221aa   , e15, f15, g15, h15, e15f, f15f, g15f, h15f, 
     *         LVS100aaai, e16, f16, g16, h16, e16f, f16f, g16f, h16f, 
c
     *         LT202aa   , e17, f17, g17, h17, e17f, f17f, g17f, h17f, 
     *         LVS101aaai, e18, f18, g18, h18, e18f, f18f, g18f, h18f, 
c
     *         Sindex, 
     *         scr13, scr14, scr15, 
     *         scr16, scr17, scr18,  
     *         nsend)
c 
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'saved_data.h'
c 
      integer e13, f13, g13, h13, e13f, f13f, g13f, h13f 
      integer e14, f14, g14, h14, e14f, f14f, g14f, h14f 
c
      integer e15, f15, g15, h15, e15f, f15f, g15f, h15f 
      integer e16, f16, g16, h16, e16f, f16f, g16f, h16f 
c
      integer e17, f17, g17, h17, e17f, f17f, g17f, h17f 
      integer e18, f18, g18, h18, e18f, f18f, g18f, h18f 
c
      integer a, b, c, i, j, k, a3, a3_s, a3_f  
      integer a1, b1, c1, i1, j1, k1  
      integer a2, b2, c2, i2, j2, k2  
      integer aa, aa1, aa2, ii, ii1, ii2 
      integer aa_range, ii_range, aa1_range, ii1_range, ii2_range, 
     *        aa2_range, a3_range   
      integer Top1(mx_array_index), Top2(mx_array_index) 
      integer ycount, ytemp, nsend,
     *        Sindex(18,mx_array_index)  
      integer y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, 
     *        y13, y14, y15, y16, y17, y18 
      integer y0, z0, z1, z2, z3, z4, z5, z6 
c
      double precision x(*), etemp, etemp1                      
c
      double precision LT102aa   ( e13:e13f,f13:f13f,g13:g13f,h13:h13f) 
      double precision LVS021aaai( e14:e14f,f14:f14f,g14:g14f,h14:h14f) 
c
      double precision LT221aa   ( e15:e15f,f15:f15f,g15:g15f,h15:h15f) 
      double precision LVS100aaai( e16:e16f,f16:f16f,g16:g16f,h16:h16f) 
c
      double precision LT202aa   ( e17:e17f,f17:f17f,g17:g17f,h17:h17f) 
      double precision LVS101aaai( e18:e18f,f18:f18f,g18:g18f,h18:h18f) 
c 
      double precision scr13(*)  
      double precision scr14(*)  
      double precision scr15(*)  
      double precision scr16(*)  
      double precision scr17(*)  
      double precision scr18(*)  
c
c ------------------------------------------------------------------------- 
c     Set output array's ranges 
c ------------------------------------------------------------------------- 
c 
      a1 = Top1(1) 
      a2 = Top2(1) 
c 
      i1 = Top1(2) 
      i2 = Top2(2) 
c 
      b1 = Top1(3) 
      b2 = Top2(3) 
c 
      j1 = Top1(4) 
      j2 = Top2(4) 
c 
      c1 = Top1(5) 
      c2 = Top2(5) 
c 
      k1 = Top1(6) 
      k2 = Top2(6) 
c
c ------------------------------------------------------------------------- 
c     Determine the contracted index ranges 
c ------------------------------------------------------------------------- 
c
      do i = 1, 4 
      do j = 1, 4 
         if (Sindex(13,i) .eq. Sindex(14,j)) then 
            if (i .eq. 1) then 
               a3_s = e13  
               a3_f = e13f  
               go to 100 
            endif 
            if (i .eq. 2) then 
               a3_s = f13  
               a3_f = f13f  
               go to 100 
            endif 
            if (i .eq. 3) then 
               a3_s = g13  
               a3_f = g13f  
               go to 100 
            endif 
            if (i .eq. 4) then 
               a3_s = h13  
               a3_f = h13f  
               go to 100 
            endif 
         endif 
      enddo 
      enddo 
100   continue 
c
c ------------------------------------------------------------------------- 
c Put arrays into scratch arrays to be used sequentially 
c ------------------------------------------------------------------------- 
c
      aa_range  = a2-a1+1 
      ii_range  = i2-i1+1 
      aa1_range = b2-b1+1 
      ii1_range = j2-j1+1 
      aa2_range = c2-c1+1 
      ii2_range = k2-k1+1 
      a3_range  = a3_f-a3_s+1 
c
      ycount = 0 
      do aa2 = c1, c2 
      do ii1 = j1, j2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y14 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii1-j1)*aa_range*a3_range  
     *         + (aa2-c1)*ii1_range*aa_range*a3_range 

             scr14(y14) = LVS021aaai(aa,a3,aa2,ii1)

            ycount = ycount + 1 
            if (y14 .ne. ycount) then 
            write(6,*) 'Y14', y14, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa1 = b1, b2 
      do ii  = i1, i2 
         do a3 = a3_s, a3_f 
            y13 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa1-b1)*ii_range*a3_range  
     *         + (ii2-k1)*aa1_range*ii_range*a3_range 

             scr13(y13) = LT102aa(aa1,ii,a3,ii2)

            ycount = ycount + 1 
            if (y13 .ne. ycount) then 
            write(6,*) 'Y13', y13, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do aa1 = b1, b2 
      do ii  = i1, i2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y16 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *          + (ii-i1)*aa_range*a3_range  
     *          + (aa1-b1)*ii_range*aa_range*a3_range 

            scr16(y16) = LVS100aaai(aa1,a3,aa,ii)

            ycount = ycount + 1 
            if (y16 .ne. ycount) then 
            write(6,*) 'Y16', y16, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa2 = c1, c2 
      do ii1 = j1, j2 
         do a3 = a3_s, a3_f 
            y15 = (a3-a3_s+1) + (ii1-j1)*a3_range   
     *         + (aa2-c1)*ii1_range*a3_range  
     *         + (ii2-k1)*aa2_range*ii1_range*a3_range 

            scr15(y15) = LT221aa(aa2,ii2,a3,ii1)

            ycount = ycount + 1 
            if (y15 .ne. ycount) then 
            write(6,*) 'Y15', y15, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii1 = j1, j2 
      do aa1 = b1, b2 
      do aa  = a1, a2 
         do a3 = a3_s, a3_f 
            y18 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (aa1-b1)*aa_range*a3_range  
     *         + (ii1-j1)*aa1_range*aa_range*a3_range 

             scr18(y18) = LVS101aaai(aa1,a3,aa,ii1)

            ycount = ycount + 1 
            if (y18 .ne. ycount) then 
            write(6,*) 'Y18', y18, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa2 = c1, c2 
      do ii  = i1, i2 
         do a3 = a3_s, a3_f 
            ycount = ycount + 1 
            y17 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa2-c1)*ii_range*a3_range  
     *         + (ii2-k1)*aa2_range*ii_range*a3_range 

c            scr17(y17) = LT202aa(aa2,ii,a3,ii2)

            if (y17 .ne. ycount) then 
            write(6,*) 'Y17', y17, ycount, a3, aa, ii, ii1  
            endif 
         enddo 
      enddo 
      enddo 
      enddo 

33    continue 
c
c ------------------------------------------------------------------------- 
c Perform calculation  
c ------------------------------------------------------------------------- 
c
      ycount = 0 
      do ii2 = k1, k2 
      do aa2 = c1, c2 
      do ii1 = j1, j2 
      do aa1 = b1, b2 
      do ii  = i1, i2 
      do aa  = a1, a2 

         etemp1 = 0 
         ycount = ycount + 1 

         do a3 = a3_s, a3_f 
c
            y14 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii1-j1)*aa_range*a3_range  
     *         + (aa2-c1)*ii1_range*aa_range*a3_range 

            y13 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa1-b1)*ii_range*a3_range  
     *         + (ii2-k1)*aa1_range*ii_range*a3_range 

            y16 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (ii-i1)*aa_range*a3_range  
     *         + (aa1-b1)*ii_range*aa_range*a3_range 

            y15 = (a3-a3_s+1) + (ii1-j1)*a3_range   
     *         + (aa2-c1)*ii1_range*a3_range  
     *         + (ii2-k1)*aa2_range*ii1_range*a3_range 

            y18 = (a3-a3_s+1) + (aa-a1)*a3_range   
     *         + (aa1-b1)*aa_range*a3_range  
     *         + (ii1-j1)*aa1_range*aa_range*a3_range 

            y17 = (a3-a3_s+1) + (ii-i1)*a3_range   
     *         + (aa2-c1)*ii_range*a3_range  
     *         + (ii2-k1)*aa2_range*ii_range*a3_range 
c
            etemp1 = etemp1 
c -------------------------------------------------------------------- 
     *        + LVS021aaai(aa,a3,aa2,ii1)*LT102aa(aa1,ii,a3,ii2)
c    *        + scr13(y13)*scr14(y14)

     *        + LVS100aaai(aa1,a3,aa,ii)*LT221aa(aa2,ii2,a3,ii1)
c    *        + scr16(y16)*scr15(y15)
c -------------------------------------------------------------------- 

     *        + LVS101aaai(aa1,a3,aa,ii1)*LT202aa(aa2,ii,a3,ii2)
c    *        + scr17(y17)*scr18(y18)
c
         enddo 

         x(ycount) = x(ycount) + etemp1  

      enddo 
      enddo 
      enddo 
      enddo 
      enddo 
      enddo 

      nsend = aa_range*ii_range*aa1_range*ii1_range*ii2_range*aa2_range    
c
c       etemp = etemp + LT001aa(aa,ii,a3,ii1)*LVS122aaai(aa1,a3,aa2,ii2)
c       etemp = etemp - LVS102aaai(aa1,a3,aa,ii2)*LT201aa(aa2,ii,a3,ii1)
c       etemp = etemp - LVS022aaai(aa,a3,aa2,ii2)*LT101aa(aa1,ii,a3,ii1)
c       etemp = etemp - LT021aa(aa,ii2,a3,ii1)*LVS120aaai(aa1,a3,aa2,ii)
c       etemp = etemp - LT002aa(aa,ii,a3,ii2)*LVS121aaai(aa1,a3,aa2,ii1)
c       etemp = etemp + LVS020aaai(aa,a3,aa2,ii)*LT121aa(aa1,ii2,a3,ii1)
c       etemp = etemp + LVS021aaai(aa,a3,aa2,ii1)*LT102aa(aa1,ii,a3,ii2)
c       etemp = etemp + LVS100aaai(aa1,a3,aa,ii)*LT221aa(aa2,ii2,a3,ii1)
c       etemp = etemp + LVS101aaai(aa1,a3,aa,ii1)*LT202aa(aa2,ii,a3,ii2)

      return
      end

