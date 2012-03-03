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
      subroutine find_non_matching_indices(x, y, array_table,
     *                             narray_table, x_unmatch, nx_unmatch,
     *                                           y_unmatch, ny_unmatch)
c------------------------------------------------------------------------
c   Find the indices of each array (x and y) that do not match with the 
c   other.
c------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer i, j, nx, ny, xval, yval
      integer x, y, x_unmatch(4), nx_unmatch, y_unmatch(4), ny_unmatch
      logical match

      nx_unmatch = 0
      ny_unmatch = 0

      nx = array_table(c_nindex,x)
      ny = array_table(c_nindex,y)
 
      do i = 1, nx
         match = .false.
         xval = array_table(c_index_array1+i-1,x)
         do j = 1, ny
            if (xval .eq. 
     *          array_table(c_index_array1+j-1,y)) then
                match = .true.   
            endif
         enddo

         if (.not. match) then
            nx_unmatch = nx_unmatch + 1
            x_unmatch(nx_unmatch) = xval
         endif
      enddo

      do j = 1, ny
         match = .false.
         yval = array_table(c_index_array1+j-1,y)
         do i = 1, nx
            if (array_table(c_index_array1+i-1,x) .eq. 
     *          yval) then
               match = .true.
            endif
         enddo

         if (.not. match) then
            ny_unmatch = ny_unmatch + 1
            y_unmatch(ny_unmatch) = yval 
         endif
      enddo

      return
      end
