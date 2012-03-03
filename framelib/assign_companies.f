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
      subroutine assign_companies(hosts, nprocs,
     *                            company,
     *                            managers_are_workers,
     *                            master_is_worker,
     *                            io_company_id)
c-------------------------------------------------------------------------
c   
c   Assigns companies and roles within the companies based on the 
c   company requirements in the "company" common block and the 
c   hardware resources in the arguments.
c
c   Arguments:
c	hosts			Host identifiers.
c	nprocs			Number of processors.
c	company			Output array containing the assigned company
c				IDs.
c       mem_assigned            Output array containing the memory assignment
c                               for each processor.
c       managers_are_workers	Flag indicating whether or not managers can
c				also be assigned as workers.
c       master_is_worker 	Flag indicating whether the master is allowed
c                               to join a company as a worker.	
c       io_company_id           ID of a special company used for I/O servers.
c                               We attempt to map one server per node if
c                               possible.
c
c-----------------------------------------------------------------------------	

      implicit none
      include 'mpif.h'
      include 'company.h'
      include 'proto_defines.h'

      integer nprocs
      integer hosts(nprocs)
      integer company(nprocs)
      integer io_company_id

      integer i, j, icompany, nm, nw, navail, memreq, ierr
      integer id_host, nmgr, nwrkr, mgr_mem, wrkr_mem, id
      integer nmapped, last_host
      integer master, pst_get_master
      logical managers_are_workers, manager_mapped
      logical master_is_worker

      do i = 1, nprocs
         company(i)      = MPI_UNDEFINED
      enddo  

      master = pst_get_master()

c---------------------------------------------------------------------------
c   Check for not enough processes.
c---------------------------------------------------------------------------

      nwrkr = 0
      do i = 1, max_company
         if (c_table(i, c_company_id) .ne. MPI_UNDEFINED) 
     *     nwrkr = nwrkr + c_table(i, c_nwrkr)
      enddo

      if (nwrkr .ne. nprocs) then
         print *,'Error: COMPANY params have requested ',nwrkr, 
     *       ' processors, but the mpirun has ',nprocs
         call mpi_abort(mpi_comm_world, 1, ierr) 
      endif

c---------------------------------------------------------------------------
c   If io_company_id .ne. 0, map the I/O company first.
c--------------------------------------------------------------------------

      if (io_company_id .ne. 0) then
         do icompany = 1, max_company
            id = c_table(icompany, c_company_id)
            if (id .eq. io_company_id) then
               nwrkr = c_table(icompany, c_nwrkr) 
               nmapped = 0

c----------------------------------------------------------------------------
c   Map "nwrkr" processes, one per host, if possible.
c----------------------------------------------------------------------------
 
   20          continue
               last_host = -1
               do i = 2, nprocs   ! never map proc 0 to io_company_id
                  if (hosts(i) .ne. last_host .and.
     *                company(i) .eq. MPI_UNDEFINED) then
                     company(i) = id     ! map the processor.
                     nmapped    = nmapped + 1
                     last_host  = hosts(i)
                     if (nmapped .eq. nwrkr) go to 50
                  endif
               enddo 
   50          continue
               
c----------------------------------------------------------------------------
c   Not enough nodes.  Continue allocating one per node until exhausted.
c----------------------------------------------------------------------------

               if (nmapped .lt. nwrkr) go to 20
            endif 
         enddo
      endif

  100 continue

c-----------------------------------------------------------------------------
c   Map the remaining companies, with as many as possible sharing the 
c   same node.
c-----------------------------------------------------------------------------

      do icompany = 1, max_company
         id = c_table(icompany, c_company_id)
         if (id .ne. io_company_id .and. 
     *       id .ne. MPI_UNDEFINED) then
            nwrkr = c_table(icompany, c_nwrkr)
            nmapped = 0
 
            do i = 1, nprocs
               if (company(i) .eq. MPI_UNDEFINED) then
                  company(i) = id    ! map the processor to this company.
                  nmapped    = nmapped + 1
                  if (nmapped .eq. nwrkr) go to 200
               endif
            enddo
  200       continue

            if (nmapped .ne. nwrkr) then
               print *,'Error: Unable to map all processors ',
     *                 ' of company ',id
               print *,'Current state of company mapping:'
               do i = 1, nprocs
                  print *,'Rank ',i,' host ',hosts(i),
     *                  ' company ',company(i)
               enddo
      
               call mpi_abort(mpi_comm_world, 1, ierr)
            endif 
         endif
      enddo

      return
      end
