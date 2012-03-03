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
      subroutine pardo_sects(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c---------------------------------------------------------------------------
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'

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

      integer ierr, array, array_type, ind
      integer type, i, j, k, nsects, vcount, ocount, vseg, oseg

c --------------------------------------------------------------------------- 
c Define the input scalar --> type of segmentation 
c --------------------------------------------------------------------------- 

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      if (array_type .ne. scalar_value) return

      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: scalar in pardo_sects, line ',
     *     current_line
         call abort_job()
      endif

      ind =  array_table(c_scalar_index, array)
      if (ind .lt. 1 .or. ind .gt. nscalar_table) then
         print *,'Scalar table index out of range in pardo_sects, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      type = scalar_table(ind)

c --------------------------------------------------------------------------- 
c Define the output scalar --> number of segmentation 
c --------------------------------------------------------------------------- 

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      if (array_type .ne. scalar_value) return

      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: scalar in pardo_sects, line ',
     *     current_line
         call abort_job()
      endif

      ind =  array_table(c_scalar_index, array)
      if (ind .lt. 1 .or. ind .gt. nscalar_table) then
         print *,'Scalar table index out of range in pardo_sects, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      oseg = scalar_table(ind)  ! The number of occupied triplets 
      vseg = nalpha_virtual/sip_mx_virt_segsize
c
c --------------------------------------------------------------------------- 
c RHF AAA triples calculation 
c --------------------------------------------------------------------------- 
c 
      if (type .eq. 1) then

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = i, vseg ! nalpha_virtual 
         do k = j, vseg ! nalpha_virtual 
            vcount = vcount + 1
         enddo
         enddo
         enddo

         nsects = my_company_size/vcount + 1
         if (oseg .lt. nsects) nsects = oseg
         if (nsects .gt. 14.0) nsects = 14.0
         if (nsects .lt. 4.0) nsects = 4.0
         scalar_table(ind) = nsects  

         return

      endif
c
c --------------------------------------------------------------------------- 
c RHF AAB triples calculation 
c --------------------------------------------------------------------------- 
c 
      if (type .eq. 2) then

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = i, vseg ! nalpha_virtual 
         do k = 1, vseg ! nalpha_virtual 
            vcount = vcount + 1
         enddo 
         enddo 
         enddo 

         nsects = my_company_size/vcount + 1
         if (oseg .lt. nsects) nsects = oseg
         if (nsects .gt. 14.0) nsects = 14.0
         if (nsects .lt. 4.0) nsects = 4.0
         scalar_table(ind) = nsects  

         return 
      
      endif
c
c --------------------------------------------------------------------------- 
c UHF aaa ring calculation 
c --------------------------------------------------------------------------- 
c 
      if (type .eq. 3) then

         vseg = nalpha_virtual/sip_mx_virt_segsize
         oseg = nalpha_occupied/sip_mx_occ_segsize

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = 1, oseg ! nalpha_occupied  
            vcount = vcount + 1
         enddo 
         enddo 

         nsects = my_company_size/vcount + 1
         if (nsects .gt. 5.0) nsects = 5.0
         scalar_table(ind) = nsects  

         return 
      
      endif
c
c ---------------------------------------------------------------------------
c UHF aaa ring calculation
c ---------------------------------------------------------------------------
c
      if (type .eq. 3) then

         vseg = nalpha_virtual/sip_mx_virt_segsize
         oseg = nalpha_occupied/sip_mx_occ_segsize

         if (vseg .lt. 1) vseg = 1
         if (oseg .lt. 1) oseg = 1

         vcount = 0
         do i = 1, vseg ! nalpha_virtual
         do j = 1, oseg ! nalpha_occupied
            vcount = vcount + 1
         enddo
         enddo

         nsects = my_company_size/vcount + 1
         if (nsects .gt. 5.0) nsects = 5.0
         scalar_table(ind) = nsects

         return

      endif

      
      
      call c_flush_stdout()
      return
      end        
         

