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
c----------------------------------------------------------------------
c
c   pst_get_master: Retrieve the master task from the pst.
c
c----------------------------------------------------------------------

      integer function pst_get_master()

      implicit none

      include 'mpif.h'
      include 'pst.h'
      include 'proto_defines.h'
      integer i, j, nprocs, me, ierr

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      do i = 1, nprocs
         if (pst(i, r_role) .eq. master_status .or.
     *       pst(i, r_role) .eq. master_worker_status) then
            pst_get_master = i - 1
            return
         endif
      enddo

      print *,'ERROR: pst_get_master cannot find the master.'
      call mpi_comm_rank(mpi_comm_world, me, ierr)
      print *,'Dump of pst on process ',me,':'
      do i = 1, nprocs
         print *,'Row ',i,' on proc ',me,' ',(pst(i,j),j=1,r_last)
      enddo 
      call mpi_abort(mpi_comm_world, 1, ierr)
      pst_get_master = -1   ! avoids HP compiler warning
      return
      end

c----------------------------------------------------------------------
c
c   pst_get_role: Retrieve the role of task "me" from the pst.
c
c----------------------------------------------------------------------

      integer function pst_get_role(me)

      implicit none
 
      include 'mpif.h'
      include 'pst.h'
      integer me, nprocs, ierr

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      if (me .ge. nprocs) then
         print *,'pst_get_role: Invalid task id: ',me
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif
      
      pst_get_role = pst(me+1, r_role)
      return
      end

c---------------------------------------------------------------------------
c   pst_get_company: Returns the company id associated with a specified
c                task.
c---------------------------------------------------------------------------

      integer function pst_get_company(iproc)
      implicit none
      include 'pst.h'
      include 'mpif.h'
      integer iproc, nprocs, ierr

      call mpi_comm_size(mpi_comm_world, nprocs, ierr) 
      if (iproc .ge. 0 .and. iproc .lt. nprocs) then
         pst_get_company = pst(iproc+1, r_company_id)
      else
         print *,'pst_get_company: Invalid processor rank = ',iproc
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif
      return
      end

c---------------------------------------------------------------------------
c   pst_get_my_company: Returns the company id associated with the running 
c                process.
c---------------------------------------------------------------------------

      integer function pst_get_my_company()
      implicit none
      include 'pst.h'
      include 'mpif.h'
      integer me, ierr

      call mpi_comm_rank(mpi_comm_world, me, ierr) 
      pst_get_my_company = pst(me+1, r_company_id)
      return
      end

c---------------------------------------------------------------------------
c   pst_get_company_comm: Returns the communicator associated with a 
c                given task's company.
c---------------------------------------------------------------------------

      integer function pst_get_company_comm(iproc)
      implicit none
      include 'pst.h'
      include 'mpif.h'
      integer iproc, nprocs, ierr

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      if (iproc .ge. 0 .and. iproc .lt. nprocs) then
         pst_get_company_comm = pst(iproc+1, r_company_comm)
      else
         print *,'pst_get_company_comm: Invalid processor rank = ',
     *            iproc
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif
      return
      end

c---------------------------------------------------------------------------
c   pst_get_company_rank: Returns the rank relative to the task's company
c                         communicator.
c---------------------------------------------------------------------------

      integer function pst_get_company_rank(iproc)
      implicit none
      include 'pst.h'
      include 'mpif.h'
      integer iproc, nprocs, ierr

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      if (iproc .ge. 0 .and. iproc .lt. nprocs) then
         pst_get_company_rank = pst(iproc+1, r_company_rank)
      else
         print *,'pst_get_company_rank: Invalid processor rank = ',
     *            iproc
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif
      return
      end


