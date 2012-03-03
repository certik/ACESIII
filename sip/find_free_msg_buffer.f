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
      integer function find_free_msg_buffer(req, nreq, itimer, ctimer)
c------------------------------------------------------------------------------
c   This function searches an array of MPI requests until it finds one that
c   is free (i. e. req(i) = MPI_REQUEST_NULL).  If it cannot find a free 
c   request, it waits for one of the requests to complete, and returns the 
c   index of the completed request.  
c-----------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'parallel_info.h'
      integer nreq
      integer req(nreq)
      integer itimer, ctimer

      integer i, ireq, ierr
      integer status(MPI_STATUS_SIZE)
   
      do i = 1, nreq
         if (req(i) .eq. MPI_REQUEST_NULL) then
            ireq = i
            go to 100
         endif
      enddo

c--------------------------------------------------------------------------
c   All buffers have outstanding requests.
c   Copy the requests into a scratch area, and test for a completion.
c--------------------------------------------------------------------------

      call update_timer(itimer)
      call timer_start(ctimer)
      call mpi_waitany(nreq, req, ireq, status, ierr)
      call update_timer(ctimer)
      call timer_start(itimer)

  100 continue
      find_free_msg_buffer = ireq
      return
      end
