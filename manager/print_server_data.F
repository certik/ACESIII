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
      subroutine print_server_data(io_company_id)
c---------------------------------------------------------------------------
c   Executed only on the master process at end of a SIAL run.  If the 
c   do_timer flag is set, we receive the accumulated server statistical
c   data from each processor in the I/O company, and print a report.
c----------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'server_stat.h'
      include 'timerz.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'interpreter.h'

      integer io_company_id
      integer lendata, niocompany
      integer ierr
      integer status(MPI_STATUS_SIZE)
      integer i, nw
      integer next

      integer pst_get_company

c---------------------------------------------------------------------------
c   Determine the number of IOCOMPANY servers.
c---------------------------------------------------------------------------

      niocompany = 0
      do i = 1, nprocs
         if (pst_get_company(i-1) .eq. io_company_id)
     *       niocompany = niocompany + 1
      enddo

      if (niocompany .eq. 0) return
      if (.not. do_timer) return

c---------------------------------------------------------------------------
c   RECV the data into the master's common block.
c--------------------------------------------------------------------------

      lendata = 12*mx_stat_keys
      call mpi_recv(sstat_tprep, lendata, MPI_DOUBLE_PRECISION, 
     *               MPI_ANY_SOURCE,
     *               sip_server_stat_data_tag, mpi_comm_world,
     *               status, ierr) 
      lendata = 6*mx_stat_keys
      call mpi_recv(sstat_nprep, lendata, MPI_INTEGER, MPI_ANY_SOURCE,
     *               sip_server_stat_data_tag, mpi_comm_world,
     *               status, ierr)

c-----------------------------------------------------------------------------
c   Print the server report.
c-----------------------------------------------------------------------------

      print *,'              SERVER PERFORMANCE STATISTICS'
      print *,'              ------ ----------- ----------' 

      print *,'Pardo      Prepares   Preparesums   Requests  Prequests'
      do i = 1, mx_stat_keys
         if (sstat_nprep(i) .gt. 0 .or. 
     *       sstat_nprepsum(i) .gt. 0 .or.
     *       sstat_nreq(i) .gt. 0 .or.
     *       sstat_npreq(i) .gt. 0) then
            print 9000,lineno(i), sstat_nprep(i),sstat_nprepsum(i),
     *                            sstat_nreq(i), sstat_npreq(i) 
         endif
      enddo

      call c_flush_stdout()

 9000 format(1x,i6,4(4x,i8)) 
      return
      end

      subroutine get_server_stats(myline, nprep, nprepsum, nreq,
     *                            npreq)
c---------------------------------------------------------------------------
c   Return server stats for a given line number.
c---------------------------------------------------------------------------
      implicit none
      include 'server_stat.h'
      integer myline, nprep, nprepsum, nreq, npreq

      integer i

      do i = 1, mx_stat_keys
         if (lineno(i) .eq. myline) then
            nprep    = sstat_nprep(i)
            nprepsum = sstat_nprepsum(i)
            nreq     = sstat_nreq(i)
            npreq    = sstat_npreq(i)
            return
         endif
      enddo

c--------------------------------------------------------------------------
c   No data for this lineno.
c--------------------------------------------------------------------------

      nprep    = 0
      nprepsum = 0
      nreq     = 0

      return
      end
