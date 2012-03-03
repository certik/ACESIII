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
      subroutine get_subrange_segments(superindex, superseg, sindex,
     *                              index_table, nindex_table,
     *                              segment_table, nsegment_table,
     *                              bseg, eseg)
c----------------------------------------------------------------------------
c   Looks up the beginning and ending segments of the sindex relative to 
c   the current segment of its superindex.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer sindex, nindex_table, nsegment_table, bseg, eseg
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)

      integer superindex, superseg, r1, r2
      integer i, r1min, r2max, r1sub, r2sub

c----------------------------------------------------------------------------
c   Look up relevant superindex info.
c----------------------------------------------------------------------------

      r1 = segment_table(c_range1,superseg)
      r2 = segment_table(c_range2,superseg)

c---------------------------------------------------------------------------
c   Search segment table.
c---------------------------------------------------------------------------

      r1min = 100000000
      r2max = 0
      bseg = 0
      eseg = 0

      do i = 1, nsegment_table
         if (segment_table(c_index,i) .eq. sindex) then
            r1sub = segment_table(c_range1,i)
            r2sub = segment_table(c_range2,i)

            if (r1sub .ge. r1 .and. r2sub .le. r2) then
               if (r1sub .le. r1min) then
                  r1min = r1sub
                  bseg  = segment_table(c_segment,i)
               endif
 
               if (r2sub .ge. r2max) then
                  r2max = r2sub
                  eseg  = segment_table(c_segment,i)
               endif
            endif
         endif
      enddo

      if (bseg .eq. 0 .or. eseg .eq. 0) then
         print *,'Error in get_subrange_segments: bseg, eseg ',
     *           bseg,eseg
         print *,'Superindex, superseg, sindex ',
     *        Superindex, superseg, sindex
         call abort_job()
      endif
 
      return
      end
 
