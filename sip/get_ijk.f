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
      subroutine get_ijk(array_table, narray_table,
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

      common /d2int_com/jatom, jx, jcenter
      common /Occ_ijk/NT, Xijk(10000,7) 
      integer jatom, jx, jcenter
      double precision flags_value

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table, sind, xind, Xijk 
      double precision scalar_table(nscalar_table), tind 
      integer*8 address_table(narray_table)
      integer*8 iarray, ievec, get_index_from_base 

      integer ierr, array, array_type, ind, nind
      integer i, nseg, bseg, eseg, start(100), end(100), maxi   
      integer Nt, wind  
      integer nindex 

      double precision x(1) 

c-----------------------------------------------------------------------
c Determine the array which will hold the 7 parameters definining
c the occupied combination.  
c------------------------------------------------------------------------
 
      array      = op(c_result_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. scalar_value) then 
         print *,'Error: First array in get_ijk must be a scalar.'
         call abort_job()
      endif

      xind =  array_table(c_scalar_index, array) 
      if (xind .lt. 1 .or. xind .gt. nscalar_table) then
         print *,'Scalar table index out of range in get_ijk, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',xind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif


c-----------------------------------------------------------------------
c Determine which combination you want   
c-----------------------------------------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array) 

      if (array_type .ne. scalar_value) then
         print *,'Error: The Second argument in get_ijk  
     *            must be a scalar.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif
 
      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in get_ijk, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

      tind = scalar_table(xind) 
      wind = scalar_table(sind)  
c     write(6,*) ' TIND XIND SIND WIND:', tind,xind, sind, wind  

      if (tind .eq. 1.0) scalar_table(xind) = Xijk(wind,1)   
      if (tind .eq. 2.0) scalar_table(xind) = Xijk(wind,2)   
      if (tind .eq. 3.0) scalar_table(xind) = Xijk(wind,3)   
      if (tind .eq. 4.0) scalar_table(xind) = Xijk(wind,4)   
      if (tind .eq. 5.0) scalar_table(xind) = Xijk(wind,5)   
      if (tind .eq. 6.0) scalar_table(xind) = Xijk(wind,6)   
      if (tind .eq. 7.0) scalar_table(xind) = Xijk(wind,7)   
c     write(6,*) ' SCALAR = ', scalar_table(xind) 

      return
      end


