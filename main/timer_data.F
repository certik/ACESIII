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
      subroutine tmanal(siofile, io_company_id)
c----------------------------------------------------------------------------
c   Gathers timer data from all servers and workers.  If the ACES_SOURCE_DIR
c   environment variable is set, the routine searches for a SIAL source 
c   file with the correct name, and produces a source-level timing analysis
c   of the gathered data.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'proto_defines.h'
      include 'proto_events.h'
      include 'timerz.h'
      include 'machine_types.h'
      include 'parallel_info.h'
#ifdef ALTIX
       include 'sheap.h'
#endif

      integer io_company_id
      character*(*) siofile

      integer i, j, nprocs_save, ierr
      integer*8 irecv
      integer*8 itimer, itstat, idesc, idesc_scr
      integer*8 itimer_avg, itimer_sd
      integer status_mpi(MPI_STATUS_SIZE), status, nints
      integer request_status(mpi_status_size,2*nprocs)
      integer timer_request(2*nprocs)
      double precision timer_stats(1) 
      integer timer_descs(1) 

      integer itdesc(1)
      equivalence (itdesc(1), tdesc(1))

      integer nw
      integer*8 c_loc64

      integer pst_get_role
#ifdef HP
      integer*8 lenshmem
#else
      integer lenshmem
#endif

#ifdef ALTIX
      pointer (iptr, timer_descs)
      pointer (dptr, timer_stats)

      iptr = ishptr
      dptr =dshptr
#else
      common timer_stats, timer_descs
#endif

c--------------------------------------------------------------------------
c   Collect the server stats, and print a summary.
c--------------------------------------------------------------------------

      call print_server_data(io_company_id)

      if (.not. do_timer) return

      nprocs_save = nprocs

c--------------------------------------------------------------------------
c   Set up the timer statistic buffer.
c--------------------------------------------------------------------------

      nw = 4*max_timers
      call mem_alloc(timer_stats, nw, bytes_per_double, itstat,
     *                .true., ierr)
      if (ierr .ne. 0) then
         print *,'Cannot alloc timer memory: itstat nw ',nw
         call abort_job()
      endif
 
      nw = max_timers
      call mem_alloc(timer_stats, nw, bytes_per_double, itimer_avg,
     *               .true., ierr)
      if (ierr .ne. 0) then
         print *,'Cannot alloc timer memory: itimer_avg nw ',nw
         call abort_job()
      endif
 
      call mem_alloc(timer_stats, nw, bytes_per_double, itimer_sd,
     *               .true., ierr)
      if (ierr .ne. 0) then
         print *,'Cannot alloc timer memory: itimer_sd nw ',nw
         call abort_job()
      endif

      nints = (max_timers*max_timer_desc_len+intsize-1)/intsize
      nw = 2*nints
      call mem_alloc(timer_descs, nw, intsize, idesc,
     *               .true., ierr)
      if (ierr .ne. 0) then
         print *,'Cannot alloc timer memory: idesc nw ',nw
         call abort_job()
      endif
      
c--------------------------------------------------------------------------
c   Initialize the sum, avg, and sd buffers.
c--------------------------------------------------------------------------

      call sum_timer_data(timer_stats(itstat), 
     *                    timers,
     *                    timer_descs(idesc),             
     *                    itdesc,
     *                    timer_stats(itimer_avg),
     *                    timer_stats(itimer_sd), max_timers,
     *                    max_timer_desc_len, .true.) 

c--------------------------------------------------------------------------
c   Post a recv for each contraction worker and integral worker in which
c   to receive the timer statistics for the workers.
c
c   These recv's will not complete until each process sends the master its
c   timer data at termination time.
c--------------------------------------------------------------------------

      do i = 2, nprocs
         status = pst_get_role(i-1)
         if (status .eq. worker_status .or.
     *       status .eq. both_status .or.
     *       status .eq. master_worker_status) then

c--------------------------------------------------------------------------
c   Receive the timer data from proc i-1.
c--------------------------------------------------------------------------

            itimer = itstat + 3*max_timers
            call mpi_recv(timer_stats(itimer), max_timers,
     *               mpi_double_precision, (i-1),
     *               timer_data_request_tag, mpi_comm_world,
     *               status_mpi, ierr)

c--------------------------------------------------------------------------
c   Receive the timer descriptor data from proc i-1.
c--------------------------------------------------------------------------

            irecv = nints + idesc   
            call mpi_recv(timer_descs(irecv), nints,
     *                  mpi_integer, (i-1),
     *                  timer_desc_request_tag, mpi_comm_world,
     *                  status_mpi, ierr)

            call sum_timer_data(timer_stats(itstat), 
     *                          timer_stats(itstat+3*max_timers),
     *                          timer_descs(idesc),             
     *                          timer_descs(idesc+nints),
     *                          timer_stats(itimer_avg),
     *                          timer_stats(itimer_sd), max_timers,
     *                          max_timer_desc_len, .false.) 
         endif
      enddo

c----------------------------------------------------------------------------
c   Print the timer data.
c----------------------------------------------------------------------------

      call timer_report(siofile, timer_descs(idesc), 
     *          timer_stats(itstat), 
     *          timer_stats(itimer_avg),
     *          timer_stats(itimer_sd),
     *          nprocs_save)

      return
      end

      subroutine sum_timer_data(sumdata, instance, descs, 
     *                          descs_instance, contrib, sumsq,
     *                          max_timers, max_timer_desc_len,
     *                          init)  
c---------------------------------------------------------------------------
c   Sums the timer data from one processor (instance) into the total sum 
c   buffer.  The number of contributing elements and the sum of squares 
c   for each timer is also accumulated.
c---------------------------------------------------------------------------
      
      implicit none
      integer i, j
      integer max_timers
      integer max_timer_desc_len
      double precision sumdata(max_timers,3), instance(max_timers)
      double precision sumsq(max_timers), contrib(max_timers)
      character*(max_timer_desc_len) descs(max_timers)
      character*(max_timer_desc_len) descs_instance(max_timers)
      logical init

      if (init) then
         do i = 1, max_timers
            descs(i) = descs_instance(i)
            if (descs(i)(1:1) .ne. ' ' .and. 
     *          descs(i)(1:1) .ne. char(0) ) then
               contrib(i) = 1.
               sumdata(i,1) = instance(i)
               sumdata(i,2) = instance(i)
               sumdata(i,3) = instance(i)
               sumsq(i) = instance(i)*instance(i)
            else
               contrib(i)   = 0.
               sumsq(i)     = 0.
               sumdata(i,1) = 0.
               sumdata(i,2) = 0.d0
               sumdata(i,3) = 1.d10 
            endif
         enddo 

         return
      endif

      do i = 1, max_timers
         if (descs_instance(i)(1:1) .ne. ' ' .and.
     *       descs_instance(i)(1:1) .ne. char(0)) then
            if (descs_instance(i) .eq. descs(i)) then
              
c------------------------------------------------------------------------
c   Sum in the data.
c------------------------------------------------------------------------
   
               if (instance(i) .ne. 0.d0) then
                  sumdata(i,1) = sumdata(i,1) + instance(i)
                  sumdata(i,2) = max(sumdata(i,2), instance(i))
                  sumdata(i,3) = min(sumdata(i,3), instance(i))
                  contrib(i)   = contrib(i) + 1.
                  sumsq(i)     = sumsq(i) + instance(i)*instance(i)
               endif   
            else
               print *,'Error: Timer desc ',descs_instance(i),
     *             ' does not match proc 0 timer_desc ',
     *              descs(i)
               print *,'i = ',i
               call abort_job()
            endif
         endif
      enddo
      return
      end
