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
      subroutine determine_array_indices(op, array_table, narray_table,
     *                                   debug)
      implicit none
      include 'interpreter.h'

      integer narray_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)

      integer x1, x2, y, opcode
      integer ny, nx1, nmatches, i, j, matches(4)
      integer x1_unmatch(4), x2_unmatch(4), nx1_unmatch, nx2_unmatch
      integer iy, ix2, x1_index

      logical debug

      opcode = op(c_opcode)
      x1     = op(c_op1_array)
      x2     = op(c_op2_array)
      y      = op(c_result_array)

      ny = array_table(c_nindex,y)
      nx1 = array_table(c_nindex,x1)
      if (debug) print *,'@@@ opcode, x1, x2, y, ny, nx1 = ',
     *             opcode, x1, x2, y, ny, nx1

      if (opcode .eq. sum_op) then
         call find_matching_indices(x1, x2, array_table, narray_table,
     *                                 matches, nmatches) 
         if (nmatches .ne. nx1) then
            print *,'Error: All indices must match for sum operation.'
            print *,'Arrays ',x1,' and ',x2,' have ',nmatches,
     *                ' matches out of ',nx1
            call abort_job()
         endif

         do i = 1, ny
            array_table(c_index_array1+i-1,y) = matches(i)
         enddo

      else if (opcode .eq. contraction_op) then
         if (debug) print *,'Before find_non_matching_indices...'
         call find_non_matching_indices(x1, x2, 
     *                                  array_table, narray_table,
     *                                  x1_unmatch, nx1_unmatch,
     *                                  x2_unmatch, nx2_unmatch)

c--------------------------------------------------------------------------
c   Step through the indices of x1, replacing each matching index with an
c   unmatched one.
c-------------------------------------------------------------------------- 

         iy = 0
         ix2 = 0
         do i = 1, nx1
            x1_index = array_table(c_index_array1+i-1, x1)
            if (debug) print *,'@@@ i = ',i,' x1_index = ',x1_index
            do j = 1, nx1_unmatch
               if (x1_index .eq. x1_unmatch(j)) then
                  iy = iy + 1
                  array_table(c_index_array1+iy-1,y) = x1_index
                  go to 100
               endif
            enddo

c--------------------------------------------------------------------------
c   This index matches an index in x2.  Replace it with an unmatching one
c   from the x2_unmatch array.
c--------------------------------------------------------------------------

            if (debug) print *,'   @@@ After j loop...'
            ix2 = ix2 + 1
            if (ix2 .gt. nx2_unmatch) go to 100
            iy  = iy + 1
            if (debug) print *,'   @@@ ix2, iy = ',ix2,iy
            array_table(c_index_array1+iy-1,y) = x2_unmatch(ix2) 
  100    continue
            if (debug) print *,'@@@ Statement 100'
         enddo
      endif

      return
      end
