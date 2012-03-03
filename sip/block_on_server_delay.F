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
      subroutine block_on_server_delay(proc)
      implicit none
      include 'mpif.h'
      include 'server_delay.h'
      include 'server_barrier_data.h'

      integer proc
#ifdef ALTIX
      integer lenshmem

      integer i, ierr
      integer status(MPI_STATUS_SIZE)
      logical done
#endif

c-------------------------------------------------------------------------
c   This subroutine is used on the Altix to block until the server delay_flag
c   on server "proc" is cleared.  This is done to allow the server time
c   to catch up with its messge processing, so it is not flooded with 
c   messages that it has no time to process.
c-------------------------------------------------------------------------

#ifdef ALTIX
      RETURN
  100 continue

c-------------------------------------------------------------------------
c   Fetch the server delay flag from server "proc".
c-------------------------------------------------------------------------

      lenshmem = 1
      call shmem_get4(delay_flag, delay_flag, lenshmem, proc)

      if (delay_flag .eq. 1) then

c-------------------------------------------------------------------------
c   Attempt to clear any leftover server_barrier requests.
c-------------------------------------------------------------------------

         do i = 1, mx_msg
            if (server_requests(i) .ne. MPI_REQUEST_NULL) then
               call mpi_test(server_requests(i), done, status, ierr)
            endif
         enddo

         call f_usleep(10000)
         go to 100
      endif 
#endif

      return
      end

      
