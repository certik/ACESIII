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
      subroutine get_index_segment(index, segment, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
c---------------------------------------------------------------------------
c   Returns the range values of a given index's segment in "val1" and "val2".
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer index, segment, nsegment_table, val1, val2
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)

      integer i, j, k, search, index_type
      integer windex

      search = index
      index_type = index_table(c_index_type, index)

      do i = 1, nsegment_table
         if (search .eq. segment_table(c_index,i) .and. 
     *       segment .eq. segment_table(c_segment,i)) then
            val1 = segment_table(c_range1,i)
            val2 = segment_table(c_range2,i)
            if ((val1 .le. 0 .or. val2 .le. 0) .or.
     *           (val1 .gt. val2)) then
               if (.not. simulator) then
                  print *,'Task ',me,' INVALID SEGMENT TABLE ENTRY'
                  print *,'  index ',search,' segment ',
     *               segment,' line ',current_line
                  do j = 1, nsegment_table
                     print *,'Entry ',i,' : ',(segment_table(k,j),
     *               k = 1, lsegment_table_entry)
                  enddo
                  call abort_job()
               endif
            endif   
            return
         endif
      enddo

c--------------------------------------------------------------------------
c   Index and segment was not found...abort job.
c--------------------------------------------------------------------------

      print *,'Task ',me,' Error: Index ',search,' and segment ',
     *        segment,' not in segment table.'
      print *,'Current instruction ',current_op,' current line ',
     *         current_line
      print *,'SEGMENT TABLE:'
      do j = 1, nsegment_table
         print *,(segment_table(i,j),i=1,lsegment_table_entry)
      enddo
      call abort_job()

      return
      end
