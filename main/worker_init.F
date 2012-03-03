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
      subroutine worker_init(istat, nworkers_needed)
c--------------------------------------------------------------------
c
c   worker_init: Handles initialization of worker tasks.
c 
c---------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'proto_defines.h'
      include 'proto_events.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'sip_tables.h'
      include 'dbugcom.h'
      include 'trace.h'
      include 'timerz.h'
#ifdef ALTIX
       include 'sheap.h'
#endif

      integer nworkers_needed
      integer master, ierr, my_company
      integer istat
      integer i
      integer*8 size, nw, indmem, i8temp

      integer pst_get_master
      integer pst_get_company  
      integer pst_get_my_company  
      integer nwmem
      integer nserver_bufs
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#endif
      integer*8 malloced_previously
      integer n, str_trimlen
      integer server_blksize
      integer niocompany

      logical heap_flag

#ifdef ALTIX
      dptr = dshptr ! set x ptr to shared memory.
      call set_memory_pointer(dshptr)
#endif
      nserver_bufs = 20

      istat = 0
      call mpi_comm_rank(mpi_comm_world, me, ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      if (dbg) print *,'Initializing worker #',me

      my_company = pst_get_my_company()
      master = pst_get_master()
      call init_timers()

      if (my_company .eq. io_company_id) then
         if (dryrun) return

c---------------------------------------------------------------------------
c   Determine how much memory remains and allocate it as temp_memory.
c---------------------------------------------------------------------------

         i8temp = maxmem * 1024
         size = i8temp * 1024
         call mem_alloc_query(malloced_previously)
         size = size - malloced_previously
         nwmem = size / bytes_per_double
         call mem_alloc(x, nwmem, bytes_per_double, indmem, .true.,
     *                     ierr)
         if (ierr .ne. 0) then
            print *,'Task ',me,' mem_alloc error in worker_init ',
     *           nwmem
            call abort_job()
         endif

         if (dbg) print *,'Actual size of main mem_alloc is ',nwmem,
     *              ' words. malloced_previously = ',
     *              malloced_previously,' words.'

c--------------------------------------------------------------------------
c   Initialize the server.
c--------------------------------------------------------------------------

         n = str_trimlen(local_path)
         server_blksize = max(sip_mx_ao_segsize,
     *                        sip_mx_occ_segsize,
     *                        sip_mx_virt_segsize)**4
         
         niocompany = 0
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) 
     *        niocompany = niocompany + 1
         enddo

         call sip_server_init(local_path(1:n), x(indmem), nwmem,
     *            nserver_bufs, server_blksize, niocompany,
     *            do_timer)
      else 
        call init_prefetch_context()   ! initializes data struct for prefetching
 
c---------------------------------------------------------------------------
c   Initialize integral computation on both servers and SIAL workers.
c---------------------------------------------------------------------------

         call sip_integral_init()

c---------------------------------------------------------------------------
c   Call the proper initialization code for contraction company.
c---------------------------------------------------------------------------

         call sip_fmain_init(istat, nworkers_needed)
      endif

      return
      end

