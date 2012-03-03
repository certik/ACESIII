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
      subroutine wait_on_block(array, block, blkndx, type, request,
     *             instruction_timer, comm_timer)
c--------------------------------------------------------------------------
c   Waits on a block that is engaged in a communication operation.
c   Handles the steps needed to properly account for the block_wait time.
c   Updates necessary block flags.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'server_monitor.h'
      include 'trace.h'
      INCLUDE 'timerz.h'

      integer array, block, type, request, instruction_timer,
     *        comm_timer, blkndx
      integer ierr
      integer status(MPI_STATUS_SIZE)
      logical flag

c---------------------------------------------------------------------------
c   Wait and measure the wait time correctly.
c---------------------------------------------------------------------------

      if (request .eq. MPI_REQUEST_NULL) return

      call update_timer(instruction_timer)
      call timer_start(comm_timer)
      call timer_start(pardo_block_wait_timer)

  100 continue
      call mpi_test(request, flag, status, ierr)
      if (.not. flag) then
         call exec_thread_server(0)
         go to 100
      endif 
      
      call update_timer(comm_timer)
      call update_timer(pardo_block_wait_timer)
      call timer_start(instruction_timer)

c--------------------------------------------------------------------------
c   Set flags properly on the block we just waited for.
c--------------------------------------------------------------------------

      if (type .ne. dummy_array_type) then
         call set_block_request(array, block, blkndx,
     *                     mpi_request_null)
         call clear_block_request_outstanding_flag(array,block,
     *                        blkndx)
         if (type .eq. served_array .and. server_monitor_on) 
     *          call server_monitor_write_log(blkndx, 'r',
     *                       status(mpi_source))

         call blkmgr_remove_block_from_list(comm_list_head, 
     *            comm_list_tail, blkndx, c_comm_list_ptr)
      endif

      return
      end
