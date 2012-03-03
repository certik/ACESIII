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
      subroutine f_acquire_pardo_lock(my_lock)
      implicit none
      include 'parallel_info.h'
     
      integer my_lock
       
      call acquire_pardo_lock(my_lock)   
      return
      end

      subroutine f_release_pardo_lock(my_lock)
      implicit none
      include 'parallel_info.h'

      integer my_lock

      call release_pardo_lock(my_lock)
      return
      end

      subroutine pardo_loadb_get_next_batch(comm, iop, index_table,
     *                        next_batch, last_batch)
      implicit none
      include 'mpif.h'
      include 'interpreter.h' 
      include 'trace.h'
      include 'parallel_info.h'
      include 'machine_types.h'
      include 'sip_tables.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
 
      integer comm, iop, next_batch, last_batch
      integer msg(4)
      integer i, ierr
      integer pardo_line
      integer request, request2
      integer status(MPI_STATUS_SIZE)
      integer*8 addr, indaddr, get_index_from_base
      integer optable_entry(loptable_entry)
      integer index_table(lindex_table_entry, *)
      integer my_batch_zero, my_last_batch_zero
      integer loop_init
      integer my_lock, batch_data(2)
      integer pardo_master, get_pardo_master
      integer cluster_size, get_pardo_cluster_size
      integer nclusters

      logical flag1, flag2

      integer icom(1)
#ifdef ALTIX
      pointer (iptr, icom)
#else
      common icom
#endif

c---------------------------------------------------------------------------
c   Send a job request to rank the pardo_master.
c---------------------------------------------------------------------------

      addr = optable_base + loptable_entry *(iop-1)*intsize
#ifdef ALTIX
      iptr = ishptr
#endif
      indaddr = get_index_from_base(addr, icom, 1)
      pardo_line = icom(indaddr+c_lineno-1)
      my_lock    = icom(indaddr+c_pardo_lock_index-1)

      msg(1) = pardo_job_request_tag
      msg(2) = pardo_line
      msg(3) = iop
      msg(4) = my_lock
      pardo_master = get_pardo_master()
      call mpi_isend(msg, 4, mpi_integer, pardo_master,
     *                     pardo_job_tag,
     *                     comm, request, ierr)

c--------------------------------------------------------------------------
c   Post a recv for the response from the thread server on proc 0.
c--------------------------------------------------------------------------

      call mpi_irecv(batch_data, 2, mpi_integer, pardo_master,
     *               pardo_job_request_tag, comm, request2,
     *               ierr)

c--------------------------------------------------------------------------
c   Clear the request on the send and recv.
c--------------------------------------------------------------------------

  100 continue
      call mpi_test(request,  flag1, status, ierr)
      if (.not. flag1) then
         call exec_thread_server(0)
         go to 100
      endif

  200 continue
      call mpi_test(request2, flag2, status, ierr)
      if (.not. flag2) then
         call exec_thread_server(0)
         go to 200
      endif

      next_batch = batch_data(1) 
      last_batch = batch_data(2)

      return
      end
