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
      subroutine worker_termination
c------------------------------------------------------------------------
c   worker_termination: Hook to give a worker a chance to finish up.
c------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'proto_events.h'
      include 'int_gen_parms.h'
      include 'timerz.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'saved_data.h'
      include 'dbugcom.h'
      include 'trace.h'

      integer master, i, dummy, ierr
      integer pst_get_my_company
      integer my_company
      integer pst_get_master
      integer nints
      integer treq1, treq2
      integer status(mpi_status_size)

      save
 
      if (.not. terminated_worker_termination) then
         terminated_worker_termination = .true.

         my_company = pst_get_my_company()
         master = pst_get_master()
         if (my_company .eq. io_company_id) then
            call prt_time('SIP server has terminated')
         else
            call prt_time('SIAL worker has terminated...')
         endif

c--------------------------------------------------------------------------
c   Send the timer data to the master.
c--------------------------------------------------------------------------

         if (me .ne. master .and. .not. dryrun .and. do_timer) then 
            call mpi_isend(timers, max_timers, mpi_double_precision,
     *              master, timer_data_request_tag, mpi_comm_world, 
     *              treq1, ierr)

            nints = (max_timers*max_timer_desc_len + intsize-1) / 
     *                intsize
            call mpi_isend(tdesc, nints, mpi_integer,
     *              master, timer_desc_request_tag, mpi_comm_world, 
     *              treq2, ierr)

            call prt_time('Waiting on timer data requests...')
            call mpi_wait(treq1, status, ierr)
            call mpi_wait(treq2, status, ierr)
            call prt_time('Completed timer data requests')
         endif
      endif

      return
      end
