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
      subroutine set_ijk_aaa(array_table, narray_table,
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
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer ierr, array, array_type, ind, nind
      integer i, nseg, bseg, eseg, start(100), end(100), maxi   

      integer Nt, np(10), ns(100,10), Xijk
      integer e, j, k, ii, jj, is, ie, js, je   

      array = op(c_result_array)
      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: Invalid array in set_ijk, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif

      nind = array_table(c_nindex, array)
      if (nind .ne. 1) then
         print *,'Error: set_ijk requires a 1-index array.'
         call abort_job()
      endif

c-----------------------------------------------------------------------
c------------------------------------------------------------------------

c     maxi  = 8 ! HARD CODED VFL 
      maxi = sip_sub_occ_segsize
      write(6,*) ' MAXI IN SET_IJK_AAA = ', maxi
 
c-----------------------------------------------------------------------
c Determine the index 
c------------------------------------------------------------------------
 
      ind   = array_table(c_index_array1,array)

      nseg = index_table(c_nsegments,ind) 
      bseg = index_table(c_bseg,ind) 
      eseg = index_table(c_eseg,ind) 

c     write(6,*) ' ------------------------------- '  
c     write(6,*) ' IND :', ind, nseg 
 
c-----------------------------------------------------------------------
c Determine the index ranges  
c np(i) = number of partitions of segment i 
c ns(i,k) =  
c------------------------------------------------------------------------

      do i = 1, nseg 
         call get_index_segment(ind, i, segment_table,
     *                     nsegment_table, index_table,
     *                     nindex_table, start(i), end(i)) 
      enddo 
c     write(6,*) '   ', ((start(i),end(i)), i=1,nseg) 

      do i = 1, nseg 
         np(i) = (end(i) - start(i) + 1)/maxi 
         if (maxi*np(i) .lt. end(i) - start(i) + 1)
     *       np(i) = np(i) + 1 
      enddo 

c     write(6,*) ' N(I):', (np(i), i=1,nseg) 

      do i = 1, nseg 
         e = 0 
         do k = 1, np(i) - 1 
            ns(i,k) = (end(i) - start(i) + 1)/np(i) 
            e = e + ns(i,k) 
         enddo 
         ns(i,np(i)) = (end(i) - start(i) + 1) - e 
      enddo 

      Nt = 0 
      ie = 0 
      do i = 1, nseg 
      do ii = 1, np(i) 

         is = ie + 1 
         ie = is + ns(i,ii) - 1 

          je = 0 
          do j = 1, nseg 
          do jj = 1, np(j) 

             js = je + 1 
             je = js + ns(j,jj) - 1 

          if (i .le. j) then 
             do k = 1, nseg 
             if (j .le. k) then 

                Nt = Nt + 1 
                Xijk(Nt,1) = i 
                Xijk(Nt,2) = is  
                Xijk(Nt,3) = ie  
                Xijk(Nt,4) = j 
                Xijk(Nt,5) = js 
                Xijk(Nt,6) = je
                Xijk(Nt,7) = k 

             endif 
             enddo 
           endif 
     
           enddo 
           enddo 

       enddo 
       enddo 

       do i = 1, 7 
       Xijk(Nt+1,i) = -1 
       enddo 

c      if (me .eq. 0) write(6,*) ' Total number of combinations:', Nt 
c      do k = 1, Nt 
c      write(6,*) '  ', (Xijk(k,i), i=1, 7) 
c      enddo 

c-----------------------------------------------------------------------
c Determine the number 
c-----------------------------------------------------------------------

      return
      end
