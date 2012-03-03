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
      subroutine set_effective_indices(array_table_entry, indices)
c--------------------------------------------------------------------------
c   Sets the index fields of a local array to specified quantities.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      integer array_table_entry(larray_table_entry)
      integer indices(*)
      integer i

      do i = 1, mx_array_index
         if (indices(i) .ne. 0) 
     *           array_table_entry(c_index_array1+i-1) = indices(i)
      enddo
      return
      end
