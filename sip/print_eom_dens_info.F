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
      subroutine print_eom_dens_info
     +                     (array_table, 
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
      integer sind, ind(mx_array_index)
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
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_result_array)
      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      if (array .eq. 0) then
         print *,'Error: The first argument in return_sval  
     *            must be an array.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      if (nindex .ne. 2) then 
         write(6,*) ' The first array in eomroot_print must have two 
     *                indices.'   
         write(6,*) ' Number given =', nindex  
         call abort_job() 
      endif 

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

c---------------------------------------------------------------------------
c   Look up its address.
c---------------------------------------------------------------------------

      block  = find_current_block(array, array_table(1,array),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,array) 
      iarray = get_block_index(array, block, stack, blkndx, x, .true.)  

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

c----------------------------------------------------------------------------
c   Locate the scalar.  
c---------------------------------------------------------------------------

      evec_array  = op(c_op1_array)
      evec_type   = array_table(c_array_type, evec_array)
      nindex_evec = array_table(c_nindex, evec_array)

      if (evec_type .ne. scalar_value) return

      if (evec_array .lt. 1 .or. evec_array .gt. narray_table) then
         print *,'Error: Invalid scalar in return_sval, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif 

      sind =  array_table(c_scalar_index, evec_array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in return_sval, ',
     *           'line ',current_line
         print *,'Index for array ',evec_array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      call doprint_eom_dens_info(x(iarray),scalar_table(sind),
     *             val1(1),val2(1),val1(2),val2(2))

      return 
      end 

      subroutine doprint_eom_dens_info(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   

      if (me .eq. 0)then
         if (sval .eq. 1) then 
            write(6,*) ' ' 
            write(6,99)
            write(6,*) ' '
            write(6,98) 'Relevant data for this root '  
            write(6,99)
            write(6,*) ' '
            do a = a1, a2
               stemp=array1(a,a)
               write(6,100) stemp, stemp*27.2113957
            enddo
         endif 

         if (sval .eq. 2) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,97) stemp
            enddo
         endif

         if (sval .eq. 3) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                              X-DIRECTION '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Right",stemp
            enddo
         endif

         if (sval .eq. 4) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Left ",stemp
            enddo
         endif

         if (sval .eq. 5) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,102) stemp
            enddo
         endif

         if (sval .eq. 6) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,103) stemp
            enddo
         endif

         if (sval .eq. 7) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                              Y-DIRECTION '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Right",stemp
            enddo
         endif

         if (sval .eq. 8) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Left ",stemp
            enddo
         endif

         if (sval .eq. 9) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,102) stemp
            enddo
         endif

         if (sval .eq. 10) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,103) stemp
            enddo
         endif

         if (sval .eq. 11) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                              Z-DIRECTION '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Right",stemp
            enddo
         endif

         if (sval .eq. 12) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,101) "Left ",stemp
            enddo
         endif

         if (sval .eq. 13) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,102) stemp
            enddo
         endif

         if (sval .eq. 14) then 
            do a = a1, a2
               stemp=array1(a,a)
               write(6,103) stemp
            enddo
         endif 

         if (sval .eq. 15) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               XX-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 16) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 17) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               YY-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 18) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 19) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               ZZ-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 20) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 21) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               XY-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 22) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 23) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               XZ-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 24) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
         endif

         if (sval .eq. 25) then
            write(6,*) ' '
            write(6,*) '                      -------------------------'
            write(6,*) '                               YZ-PLANE   '
            write(6,*) '                      -------------------------'
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Electric Polarizability:           ",stemp
            enddo
         endif

         if (sval .eq. 26) then
            do a = a1, a2
               stemp=array1(a,a)
               write(6,104) "Cumulative Electric Polarizability:",stemp
            enddo
            write(6,*) ' '
            write(6,99)
            write(6,99)
            write(6,*) ' '
         endif

      endif

 97   format(5x,"Oscillator Strength:",3x,F12.8)
 98   format(23x,A27)
 99   format(2x,66("-"))
100   format(5x,"Excitation Energy:",5x,F12.8," a.u.",5x,F12.8," eV")
101   format(5x, A5," Transition Moment:",23x,F12.8)
102   format(5x,"Dipole Strength:",31x,F12.8)
103   format(5x,"Oscillator Strengh:",28x,F12.8)
104   format(5x,A35,12x,F12.8)

      return
      end

