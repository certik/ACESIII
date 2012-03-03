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
      subroutine send_master_barrier_msg()
c-----------------------------------------------------------------------------
c     Sends a "msg_master_barrier" msg to each other worker in the master's 
c     company, and waits for completion of all the messages.  This messaging
c     is done by the master, after he has received a message from each 
c     worker indicating that it has reached a sip_barrier or server_barrier.
c-----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      integer statuses(MPI_STATUS_SIZE,my_company_size)      
      integer requests(my_company_size)
      integer msg(5)
      integer i
      integer ierr
      integer company_comm, pst_get_company_comm

      if (dbg) then
         print *,'MASTER: send_master_barrier_msg() was entered.'
         call c_flush_stdout()
      endif

      if (my_company_size .eq. 1) return

      company_comm = pst_get_company_comm(me)

c-----------------------------------------------------------------------------
c   Send out the messages to each worker in my company.
c-----------------------------------------------------------------------------

      msg(1) = -8   ! signal is -8, send with tag 3456.
      do i = 1, my_company_size
         if (i-1 .ne. my_company_rank) then
            call mpi_isend(msg, 1, MPI_INTEGER, i-1,
     *                     3456, company_comm, requests(i), ierr)
         else
            requests(i) = MPI_REQUEST_NULL
         endif
      enddo

c-----------------------------------------------------------------------------
c   Wait for completion of the messages.
c-----------------------------------------------------------------------------

      if (dbg) print *,'Task ',me,' Master sent master_barrier_msg ',
     *   ' to each worker.'
      call mpi_waitall(my_company_size, requests, statuses, ierr)
      if (dbg) print *,'Task ',me,
     *   ' After mpi_waitall in send_master_barrier_msg, line ',
     *                     current_line
      return
      end
