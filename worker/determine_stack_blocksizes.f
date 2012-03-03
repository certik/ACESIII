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
      subroutine determine_stack_blocksizes(index_table, nindex_table,
     *               segment_table, nsegment_table, boccval, eoccval,
     *               baoccval, eaoccval, bboccval, eboccval,
     *               stack_blocksizes, nstacks)
c---------------------------------------------------------------------------
c   Calculates the blocksize for each stack used.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'dbugcom.h'
      include 'sial_config_params.h'

      integer nindex_table, nsegment_table, nstacks
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer stack_blocksizes(nstacks)

      integer i, j, allmax, occmax, virtmax, aomax
      integer baoccval, eaoccval
      integer bboccval, eboccval
      integer boccval, eoccval
      integer iseg, range, index, index_type
      integer minblk, imin, temp

c----------------------------------------------------------------------------
c   Find the maximum segment size of each tye of index.
c----------------------------------------------------------------------------

      allmax = 0
      occmax = 0
      virtmax = 0
      aomax = 0

      do i = 1, nsegment_table
         index = segment_table(c_index,i)
         index_type = index_table(c_index_type,index)

         iseg   = segment_table(c_segment,i)
         range  = segment_table(c_range2,i)-segment_table(c_range1,i)+1
         allmax = max(allmax,range)
         if (index_type .eq. moaindex) then 
            if (iseg .ge. baoccval .and. iseg .le. eaoccval) then
               occmax = max(occmax,range)
            else
               virtmax = max(virtmax,range)
            endif
         else if (index_type .eq. mobindex) then
            if (iseg .ge. bboccval .and. iseg .le. eboccval) then
               occmax = max(occmax,range)
            else
               virtmax = max(virtmax,range)
            endif 
         else if (index_type .eq. moindex) then
            if (iseg .ge. boccval .and. iseg .le. eoccval) then
               occmax = max(occmax,range)
            else
               virtmax = max(virtmax,range)
            endif
         else if (index_type .eq. aoindex) then
            aomax = max(aomax, range)
         endif   
      enddo

c----------------------------------------------------------------------------
c   Normal stack allocation:
c   Stack 1: 2-dimensional arrays.
c   Stack 2: OOOV/OOOO
c   Stack 3: Interpolated between 2 & 4
c   Stack 4: OOVV
c   Stack 5: Interpolated between 4 & 6
c   Stack 6: OVVV
c   Stack 7: VVVV
c
c   With VVVI_STACK = .true.
c   Stack 1: 2-dimensional arrays.
c   Stack 2: OOOV/OOOO
c   Stack 3: Interpolated between 2 & 4
c   Stack 4: VVVI
c   Stack 5: Interpolated between 4 & 6
c   Stack 6: OOVV
c   Stack 7: Interpolated between 6 & 8 
c   Stack 8: OVVV
c   Stack 9: VVVV
c----------------------------------------------------------------------------

      virtmax = max(virtmax, aomax)

      if (dbg) print *,'allmax, occmax, virtmax = ',
     *                  allmax,occmax,virtmax
      if (vvvi_stack) then
         if (nstacks .lt. 9) then
            print *,'Error: VVVI_STACK is set in sial_config.  '
            print *,'       Requires at least 9 memory stacks.'
            call abort_job()
         endif

         stack_blocksizes(1) = allmax * allmax
         stack_blocksizes(2) = max(occmax*occmax*occmax*occmax,
     *                          occmax*occmax*occmax*virtmax)
         stack_blocksizes(4) = virtmax*virtmax*virtmax
         stack_blocksizes(6) = occmax*occmax*virtmax*virtmax
         stack_blocksizes(8) = occmax*virtmax*virtmax*virtmax
         stack_blocksizes(9) = virtmax*virtmax*virtmax*virtmax

         stack_blocksizes(3) = 
     *             (stack_blocksizes(2)+stack_blocksizes(4))/2
         stack_blocksizes(5) = 
     *             (stack_blocksizes(6)+stack_blocksizes(4))/2
         stack_blocksizes(7) = 
     *             (stack_blocksizes(6)+stack_blocksizes(8))/2
      else 
         if (nstacks .lt. 7) then
            print *,'Error: At least 7 memory stacks are required.'
            call abort_job()
         endif
         stack_blocksizes(1) = allmax * allmax
         stack_blocksizes(2) = max(occmax*occmax*occmax*occmax,
     *                          occmax*occmax*occmax*virtmax)
         stack_blocksizes(4) = occmax*occmax*virtmax*virtmax
         stack_blocksizes(6) = occmax*virtmax*virtmax*virtmax
         stack_blocksizes(7) = virtmax*virtmax*virtmax*virtmax

         stack_blocksizes(3) = 
     *             (stack_blocksizes(2)+stack_blocksizes(4))/2
         stack_blocksizes(5) = 
     *             (stack_blocksizes(6)+stack_blocksizes(4))/2
      endif

      if (dbg) print *,'stack blocksizes: ',
     *                   (stack_blocksizes(i),i=1,nstacks)

c--------------------------------------------------------------------------
c   Sort the stacks into increasing order.
c--------------------------------------------------------------------------

      do i = 1, nstacks
         minblk = stack_blocksizes(i)
         imin = i
         do j = i+1,nstacks
            if (stack_blocksizes(j) .lt. minblk) then
               minblk = stack_blocksizes(j)
               imin = j
            endif
         enddo

         if (imin .ne. i) then
            temp = stack_blocksizes(i)
            stack_blocksizes(i) = stack_blocksizes(imin)
            stack_blocksizes(imin) = temp
         endif
      enddo

      if (dbg) print *,'Sorted stack_blocksizes: ',
     *     (stack_blocksizes(i),
     *     i = 1, nstacks)

      return
      end
