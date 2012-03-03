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
      subroutine build_pst(nprocs, master_is_worker, 
     *                     managers_are_workers, maxmem,
     *                     io_company_id)

c------------------------------------------------------------------------
c   build_pst: Builds the pst (process status table).
c   
c   1. Assigns the master task.
c   2. The master task broadcasts a hostid request to each task.
c   3. The master obtains the hostid responses, and assigns each task
c      one or more roles (i. e. worker, manager, or both).
c   4. Communicators are built for each set of tasks on teh same host.
c   5. A communicator is built containing all the managers.
c   6. The pst is built from this data and distributed to each task.
c
c-----------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'proto_defines.h'
      include 'proto_events.h'
      include 'pst.h'
      include 'company.h'
      include 'dbugcom.h'

      integer me, nprocs
      integer master
      integer ihosts(nprocs), task_status(nprocs)
      integer io_company_id

      integer i, j, ierr, local_rank
      integer ihost_id
      integer local_comm
      integer nlocal_grp, my_manager
      integer next, len

      integer id, nmgr, mgr_mem, nwrkr, wrkr_mem
      integer maxmem
      integer mpi_group_world, ngroup, new_group
      integer icompany, my_role, mem
      integer nw
      integer company_comm, company_rank
      integer lenmode, nproc_per_node, nodes
      integer pst_get_company_comm
      integer sip_pst_tag
      parameter (sip_pst_tag = 9999)

      logical master_is_worker
      logical managers_are_workers

      character*(MPI_MAX_PROCESSOR_NAME ) hname
      character*(MPI_MAX_PROCESSOR_NAME ) hosts(nprocs)
      character*256 modevar

      character*15 role(5)
      data role/'Worker', 'Manager', 'Manager/Worker',
     *          'Master', 'Master/Worker'/
      master = 0
      call mpi_comm_rank(mpi_comm_world, me, ierr)

c-----------------------------------------------------------------------
c   Get the host names from each processor on the master.
c-----------------------------------------------------------------------

      hname = ' '
      call mpi_get_processor_name(hname, len, ierr)
      call mpi_gather(hname, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, 
     *                hosts, MPI_MAX_PROCESSOR_NAME, 
     *                MPI_CHARACTER, master, mpi_comm_world, 
     *                ierr)
      if (ierr .ne. 0) then
         print *,'Processor ',me,' ERROR IN MPI_GATHER: ',ierr
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

c-----------------------------------------------------------------------
c   Form a list of unique host identifiers from the processor names.
c-----------------------------------------------------------------------

      do i = 1, nprocs
         ihosts(i) = 0
      enddo

      next = 0
      do i = 1, nprocs
         if (ihosts(i) .eq. 0) then
            next = next + 1
            hname = hosts(i)

            do j = i, nprocs
               if (hname .eq. hosts(j)) then
                  ihosts(j) = next
               endif
            enddo
         endif 
      enddo

#ifdef BLUEGENE
      call c_getenv('mode'//char(0), modevar, lenmode, ierr)
      if (ierr .eq. 0) then
         if (modevar(1:3) .eq. 'smp') then
            nproc_per_node = 1
         else if (modevar(1:4) .eq. 'dual') then
            nproc_per_node = 2
         else if (modevar(1:2) .eq. 'vn') then
            nproc_per_node = 4
         endif 
      else
         print *,'Error: mode environment variable is required ',
     *           'on BLUE GENE.'
         call abort_job()
      endif 

      next = 0
      nodes = nprocs / nproc_per_node
      if (nprocs .gt. nodes * nproc_per_node) nodes = nodes + 1
      DO I = 1, nodes
         do j = 1, nproc_per_node
            next = next + 1 
            if (next .le. nprocs) ihosts(next) = i 
         enddo 
      ENDDO
#endif

c-----------------------------------------------------------------------
c   Assign each task its status.
c-----------------------------------------------------------------------

      if (me .eq. master) then
         call proto_assign_status(master, ihosts, nprocs, task_status,
     *                    master_is_worker, managers_are_workers)

         if (dbg) then
            print *,'----------- Prototype setup -----------------'
            print *,'Number of tasks.................',nprocs
            print *,'Rank of master task.............',master
            print *,'Assigned status of tasks:'
            do i = 1, nprocs
	       print *,'     Host Name: ',
     *              hosts(i)(1:MPI_MAX_PROCESSOR_NAME),
     *              role(task_status(i))
            enddo
         endif
      endif

      call mpi_barrier(mpi_comm_world, ierr)
 
c-----------------------------------------------------------------------
c   Broadcast the ihosts and status arrays.
c----------------------------------------------------------------------- 

      call mpi_bcast(task_status, nprocs, MPI_INTEGER, master, 
     *               mpi_comm_world, ierr)
      if (ierr .ne. 0) then
         print *,'Broadcast of status array failed.'
         print *,'ierr ',ierr
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(ihosts, nprocs, MPI_INTEGER, master, 
     *               mpi_comm_world, ierr)
      if (ierr .ne. 0) then
         print *,'Broadcast of ihosts array failed.'
         print *,'ierr ',ierr
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

c------------------------------------------------------------------------
c   Fill in the "role"  and "status" fields of the pst.
c------------------------------------------------------------------------

      do i = 1, nprocs
         pst(i, r_role)   = task_status(i)
      enddo

c-----------------------------------------------------------------------
c   Create a local comm for each node.
c-----------------------------------------------------------------------

      nlocal_grp = 0
      ihost_id = ihosts(me+1)
      my_manager = 0
      do i = 1, nprocs
         if (ihosts(i) .eq. ihost_id) then
            nlocal_grp = nlocal_grp + 1

c------------------------------------------------------------------------
c   Find the manager of this group.
c------------------------------------------------------------------------

            if (task_status(i) .eq. manager_status .or.
     *          task_status(i) .eq. both_status) then
               if (master .ne. i-1)
     *            my_manager = nlocal_grp-1  ! rank within this grp
            endif
         endif
      enddo

c-------------------------------------------------------------------------
c   Perform the company/platoon assignments.
c-------------------------------------------------------------------------

      if (me .eq. master) then
#ifdef BLUEGENE
         call bluegene_assign_companies(ihosts, nprocs, nproc_per_node,
     *                  pst(1,r_company_id), pst(1,r_platoon_id))
#else
         call assign_companies(ihosts, nprocs, 
     *                         pst(1,r_company_id), 
     *                         managers_are_workers,
     *                         master_is_worker,
     *                         io_company_id)
#endif
      endif 

      call mpi_barrier(mpi_comm_world, ierr)
 
c-------------------------------------------------------------------------
c   Broadcast the table as completed to this point.
c-------------------------------------------------------------------------

      call mpi_bcast(pst, r_last*maxprocs, MPI_INTEGER, master,
     *               mpi_comm_world, ierr)
      if (ierr .ne. 0) then
         print *,'Broadcast of pst failed.'
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

c-------------------------------------------------------------------------
c   Form an MPI communicator for each of the companies.
c-------------------------------------------------------------------------

      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_comm_split(mpi_comm_world, pst(me+1,r_company_id),
     *                    me, pst(me+1,r_company_comm), ierr)
      company_comm = pst(me+1,r_company_comm)   ! save communicator id

      if (ierr .ne. 0) then
         print *,'MPI_comm_split on company_id failed.'
         print *,'ierr = ',ierr
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

c--------------------------------------------------------------------------
c   Gather the rank within company to the master.
c--------------------------------------------------------------------------

      if (pst(me+1,r_company_comm) .eq. MPI_UNDEFINED) then
         company_rank = MPI_UNDEFINED
      else
        call mpi_comm_rank(pst(me+1,r_company_comm), 
     *                   company_rank, ierr)
      endif
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_gather(company_rank, 1, MPI_INTEGER, 
     *                pst(1,r_company_rank), 1, MPI_INTEGER,
     *                master, mpi_comm_world, ierr)
      
c         if (me .eq. 0) then
c            print *,'NODE-TO-RANK MAPPING'
c            do i = 1, nprocs
c               print *,'Rank ',i-1,' node ',ihosts(i),' company ',
c     *            pst(i,r_company_id),' company rank ',
c     *             pst(i,r_company_rank) 
c            enddo
c         endif

      if (ierr .ne. 0) then
         print *,'Gather on company_rank failed.'
         print *,'ierr ',ierr
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

      call mpi_barrier(mpi_comm_world, ierr)
      if (me .eq. 0) print *,'Gather on company_rank succeeded.'

      if (dbg .and. me .eq. master) then
         print *,'------------ Initial PST table ------------'
         do i = 1, nprocs
            print *,(pst(i,j),j=1,r_last)
         enddo
      endif
      
c-------------------------------------------------------------------------
c   Broadcast the completed table.
c-------------------------------------------------------------------------

      call mpi_bcast(pst, r_last*maxprocs, MPI_INTEGER, master, 
     *               mpi_comm_world, ierr)
      if (ierr .ne. 0) then
         print *,'Broadcast of pst failed.'
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

      call mpi_barrier(mpi_comm_world, ierr)
      pst(me+1,r_company_comm) = company_comm  ! store comm id in table. 
 
      return
      end
