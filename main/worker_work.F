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
      subroutine worker_work()
c---------------------------------------------------------------------------
c   "Work" routine for the worker task.
c---------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'proto_events.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'trace.h'
      include 'saved_data.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer pst_get_my_company
      integer pst_get_master
      integer my_company
      integer ierr, master
      integer dummy
      integer dryrun_flag
      integer dbgval

      double precision t1
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer msg1, msg2
      save msg1, msg2

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      my_company = pst_get_my_company()
c#ifdef HP
      call clear_saved_requests()
c#endif

      if (dryrun) return

#ifdef ALTIX
      dptr = dshptr
#endif

c-------------------------------------------------------------------------
c   Call worker subroutines for the particular company.
c-------------------------------------------------------------------------

      if (my_company .eq. io_company_id) then
         call sip_server_work()   ! enter server processing
      else if (my_company .ne. MPI_UNDEFINED) then

         dryrun_flag = 0

         if (dbg) then
            dbgval = 1
         else
            dbgval = 0
         endif

         call sumz_work(dryrun_flag, nworkthread, dbgval, x(ithread))
         if (dbg) print *,'Task ',me,' Exiting sumz_work...'
      endif
      return
      end
