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
      subroutine find_matching_indices(x, y, array_table, narray_table,
     *                                 matches, nmatches)
c---------------------------------------------------------------------------
c   Returns matching indices of arrays x and y in array "matches".
c---------------------------------------------------------------------------
 
      implicit none
      include 'interpreter.h'
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer x, y, matches(4)
      integer i, j, nmatches, nx, ny, xval

      nx = array_table(c_nindex,x)    ! number of indices.
      ny = array_table(c_nindex,y)

      nmatches = 0
      do i = 1, nx
         xval = array_table(c_index_array1+i-1,x)
         do j = 1, ny
            if (xval .eq. array_table(c_index_array1+j-1,y)) then
               nmatches = nmatches + 1
               matches(nmatches) = xval
            endif
         enddo 
      enddo
      
      return
      end
