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
      subroutine copy_ba(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c
c   The elements of Array1 are put into Array2 for overlaping indeces. 
c
c   The instruction is executed as follows:
c   
c   execute copy_ab Array1 Array2
c
c   where Array1 is the full input matrix. 
c         Array2 is the same ar Array1 where indeces overlap.  
c
c   Both Array1 and Array2 must be "static" arrays.
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
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer*8 iarray, ievec, get_index_from_base
      integer i, j, na, nb  
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer inda(mx_array_index)
      integer indb(mx_array_index)
      integer junk

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
c   Locate the data for the input matrix.
c---------------------------------------------------------------------------

      array      = op(c_result_array)
      evec_array = op(c_op1_array)

      if (array .eq. 0 .or. evec_array .eq. 0) then
         print *,'Error: remove_diagonal routine requires 2 
     *            array arguments.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      array_type = array_table(c_array_type, array)
      evec_type  = array_table(c_array_type, evec_array)

      if (array_type .ne. static_array .or. 
     *    evec_type .ne. static_array) then
         print *,'Error: Both arrays in remove_diagonal routine 
     *            must be static.'
         call abort_job()
      endif

      nindex      = array_table(c_nindex, array)
      nindex_evec = array_table(c_nindex, evec_array)

      if (nindex .ne. nindex_evec) then
         print *,'Error: Both arrays in remove_diagonal must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex,' indices.'
         print *,'Second array has ',nindex_evec,' indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         nb = array_table(c_index_array1+i-1,evec_array)

         inda(i) = index_table(c_index_size, na)  ! pick up length of index
         indb(i) = index_table(c_index_size, nb)  ! pick up length of index
 
      enddo

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray = get_index_from_base(address_table(array), x, 2)
      ievec  = get_index_from_base(address_table(evec_array), x, 2)

c     na = index_table(c_index_size, ind(1))  ! pick up length of index
c     nb = index_table(c_index_size, ind(1))  ! pick up length of index

      if (nindex .eq. 2) call copy2ba(x(iarray),x(ievec),nindex,
     *                              inda(1),inda(2),indb(1),indb(2)) 

      return 
      end 

      subroutine copy2ba(array1,array2,nindex,na1,na2,nb1,nb2) 
      implicit none

      include 'int_gen_parms.h'

c----------------------------------------------------------------------------
c
c The array1 is assumed to be of type beta-beta. The elements 
c of the array1 are put into array2. 
c 
c----------------------------------------------------------------------------

      integer nindex, na1, na2, nb1, nb2  
      double precision array1(na1,na2)
      double precision array2(nb1,nb2)
      integer a, b, n1, n2, noff   

      n1 = min0(na1,nb1) 
      n2 = min0(na2,nb2) 

      write(6,*) ' RANGE OF INDICES IN COPY2_BA ' 
      write(6,*) ' Array1 :', na1, na2   
      write(6,*) ' Array2 :', nb1, nb2   

      noff = nalpha_occupied - nbeta_occupied 

      write(6,*) ' NOFF = ', noff 

      if (na1 .eq. nalpha_occupied) then  

         noff = nbeta_occupied  
         write(6,*) ' NOFF = ', noff 

         do b = 1, n2   
            do a = 1, nalpha_occupied - nbeta_occupied    
               array2(a,b) = array2(a,b) + array1(a+noff,b)
               write(6,*) ' COPY2_BA :', a, b, array1(a+noff,b), 
     *                                      array2(a,b) 
            enddo
         enddo

         go to 100 

      endif 

      do b = 1, n2  
         do a = 1, noff 
            array2(a,b) = 0.0d0 
         enddo 
         do a = 1, n1  
            array2(a+noff,b) = array2(a+noff,b) + array1(a,b)
            write(6,*) ' COPY2_BA :', a+noff, b, array1(a,b), 
     *                                        array2(a,b) 
         enddo
      enddo

100   continue 

      return
      end

