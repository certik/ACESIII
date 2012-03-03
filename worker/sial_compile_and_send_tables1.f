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
      subroutine sial_compile_and_send_tables1()
c----------------------------------------------------------------------------
c   This subroutine does the following:
c   Master:   
c      1. Compiles the SIAL programs used for each company.
c      2. Sends the descriptions of the compiled tables to each processor 
c         in the company.
c      3. Sends the relevant tables to each processor in the company.
c
c   Worker:
c      1. Receives the table description from the master.
c      2. Creates the tables.
c      3. Receives the table data from the master.
c----------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'company.h'
      include 'int_gen_parms.h'
      include 'proto_defines.h'

      integer master
      integer pst_get_master
      integer me, nprocs, ierr

      integer company, target_company
      integer i, j, k, n
      integer interpretface, str_trimlen
      integer pst_get_company
      integer req(1000), status(mpi_status_size, 1000)
      integer wait_status(mpi_status_size)
      integer ireq, nreq

      character*80 sial_file
      call mpi_comm_rank(mpi_comm_world, me, ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      master = pst_get_master()

      call load_pre_defined_routines()

      if (master .eq. me) then
         do i = 1, max_company
            if (c_table(i,c_company_id) .ne. MPI_UNDEFINED) then 

c---------------------------------------------------------------------------
c   Compile the SIAL code (or read an object file) for this company.
c---------------------------------------------------------------------------

              sial_file = company_sial_prog(i)
              if (sial_file(1:1) .eq. ' ') go to 100 ! no program
              target_company   = c_table(i,c_company_id)

              n = str_trimlen(sial_file)
              if (sial_file(n-3:n) .eq. '.sio') then
                 call read_tables(sial_file(1:n) // char(0))
              else if (sial_file(n-4:n) .eq. '.sial') then
                  print *,'Error: Cannot use .sial files.'
                  call abort_job()
c                 ierr = interpretface(sial_file(1:n) // char(0) )
c                 if (ierr .ne. 0) then
c                    print *,'Error: SIAL compile errors.'
c                    call abort_job()
c                 endif
              else
                 print *,'Error: Invalid file format for SIAL program.'
                 print *,'File extension must be .sial or .sio'
                 print *,'SIAL_FILE = ',sial_file(1:n)
                 call abort_job()
              endif

c--------------------------------------------------------------------------
c   Distribute the tables to all processors in the company.
c--------------------------------------------------------------------------

              nreq = 0
              do j = 1, nprocs
                 company = pst_get_company(j-1)
                 if (company .eq. target_company) then
                    if (nreq .eq. 1000) then
                       do k = 1, 5   ! Need 5 empty slots per call 

c---------------------------------------------------------------------------
c   We have exhausted the amount of requests.  Wait for completion of one
c   and move req(nreq) to it's slot.
c---------------------------------------------------------------------------

                          call mpi_waitany(nreq, req, ireq, 
     *                                  wait_status, ierr)
                          if (ireq .ne. nreq) then
                             req(ireq) = req(nreq)  ! move last req down
                          endif
                          nreq = nreq - 1
                       enddo
                    endif

                    call send_tables(j-1, req, nreq)
                 endif 
              enddo

c---------------------------------------------------------------------------
c   Free the SIP tables.
c---------------------------------------------------------------------------

              call mpi_waitall(nreq, req, status, ierr)

              company = pst_get_company(me)
              if (company .eq. io_company_id)
     *             call free_sip_tables()
            endif
  100       continue
         enddo 
      else

c---------------------------------------------------------------------------
c   Receive the worker's tables from the master.
c---------------------------------------------------------------------------

         company = pst_get_company(me)
         if (company .ne. io_company_id) 
     *      call recv_tables(master) 
      endif

      return
      end
