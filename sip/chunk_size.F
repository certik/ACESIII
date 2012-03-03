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
      integer function chunk_size(next_batch_start, max_batch, np,
     *                            last_chunk)
c---------------------------------------------------------------------------
c   Determines the next chunk size for the load-balancing parallel loop.
c
c   Arguments:
c      next_batch_start - The beginning of the next chunk of work to be
c                         performed.
c      max_batch        - The total number of loop iterations to be performed.
c      np               - The number of processors executing the loop.
c----------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      integer next_batch_start, max_batch, np
      integer n, npcut
      integer last_chunk
      integer max_chunk, avg_batch
      parameter (max_chunk = 100)
      integer cutoff
      integer distribute_batch
      double precision slope, dchunk

      if (max_batch .le. np) then
         chunk_size = 1
         return
      endif

      n = max_batch - next_batch_start + 1   ! number of iterations left
      if (n .le. 0) then
         chunk_size = 1
         return
      endif

      distribute_batch = 0.9 * max_batch  ! feed out 90% of work in first pass
      avg_batch = (distribute_batch + np -1) / np
      cutoff = avg_batch * np
      if (next_batch_start .lt. cutoff) then
         chunk_size = min(avg_batch, max_chunk)
      else
            slope = real(next_batch_start-cutoff)/(max_batch-cutoff)
            dchunk = avg_batch + (1.0-avg_batch)*slope
             chunk_size = dchunk + 0.5
            chunk_size = max(min(chunk_size,max_chunk), 1)
c         endif
      endif

      chunk_size = min(chunk_size, n)
      return
      end
 
