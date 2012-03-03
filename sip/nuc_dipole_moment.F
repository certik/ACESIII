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
      subroutine nuc_dipole_moment(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)

c--------------------------------------------------------------------------
c
c The scalar1 is raised to the scalar2 power and replaced. 
c
c scalar1 = scalar1**scalar2 
c   
c--------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'
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
      integer natoms, imoment
      integer i, j, n1, n2, n3, n4  
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer sind, xind, ind(mx_array_index)
      integer seg, index(mx_array_index), val1(mx_array_index), 
     *        val2(mx_array_index), junk

      integer block, blkndx
      integer find_current_block
      integer stack

      double precision DNUC,XCOORD,YCOORD,ZCOORD
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
c   Locate the scalar that tells whether X Y or Z nuclear dipole moment
c---------------------------------------------------------------------------

      array      = op(c_result_array)
      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         print *,'Error: The argument in nn_dipole_moment
     *            must be a scalar.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
        print *,'Scalar table index out of range in nn_dipole_moment, ',
     *          'line ',current_line
        print *,'Index for array ',array,' is ',sind,' should be ',
     *          'between 1 and ',nscalar_table
        call abort_job()
      endif

c----------------------------------------------------------------------------
c   Locate the scalar that will hold the nuclear dipole moment
C
C   !!!  On input, this scalar has the number of atoms  !!!
c---------------------------------------------------------------------------

      evec_array  = op(c_op1_array)
      evec_type   = array_table(c_array_type, evec_array)
      nindex_evec = array_table(c_nindex, evec_array)

      xind =  array_table(c_scalar_index, evec_array)
      if (xind .lt. 1 .or. xind .gt. nscalar_table) then
        print *,'Scalar table index out of range in nn_dipole_moment, ',
     *          'line ',current_line
        print *,'Index for array ',evec_array,'is',xind,' should be ',
     *          'between 1 and ',nscalar_table
        call abort_job()
      endif
C
C
C   Compute the nuclear contribution to the dipole moment
C
C 
      if ( abs (scalar_table(sind) - 1) .LT. 0.0001 ) imoment = 1
      if ( abs (scalar_table(sind) - 2) .LT. 0.0001 ) imoment = 2
      if ( abs (scalar_table(sind) - 3) .LT. 0.0001 ) imoment = 3

      natoms = scalar_table(xind)
      scalar_table(xind) = 0.0D0
      DNUC = 0.0D0
      IF (imoment .EQ. 1) THEN
         DO I = 1, natoms
            XCOORD = acenter(I,1)
            DNUC = DNUC + charge(I)*XCOORD
         ENDDO ! I
      ELSE IF (imoment .EQ. 2) THEN
         DO I = 1, natoms
            YCOORD = acenter(I,2)
            DNUC = DNUC + charge(I)*YCOORD
         ENDDO ! I
      ELSE IF (imoment .EQ. 3) THEN
         DO I = 1, natoms
            ZCOORD = acenter(I,3)
            DNUC = DNUC + charge(I)*ZCOORD
         ENDDO ! I
      ELSE
         WRITE(*,*) "Wrong type of moment; should be 1, 2, or 3."
      ENDIF
      scalar_table(xind) = DNUC

      return 
      end 


