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
      subroutine sial_compile_and_send_tables()
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
      include 'sip_tables.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer create_index_table
      integer create_array_table
      integer create_optable

      integer master
      integer me, nprocs, ierr

      integer company_comm, target_company
      integer company_rank, rank
      integer i, j, n
      integer interpretface, str_trimlen
      integer pst_get_company
      integer pst_get_company_comm
      integer pst_get_company_rank
      integer msg(4)
      integer*8 ind, get_index_from_base
      integer ix(1)
      double precision x(1) 
#ifdef ALTIX
      pointer (iptr, ix)
      pointer (dptr, x)
#endif

      character*80 sial_file

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      company_comm = pst_get_company_comm(me)
      company_rank = pst_get_company_rank(me)
      target_company = pst_get_company(me)
      if (target_company .eq. io_company_id) return  ! this is a  server

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c--------------------------------------------------------------------------
c   Find the master of this company (company_rank = 0).
c--------------------------------------------------------------------------

      do i = 1, nprocs
         if (pst_get_company(i-1) .eq. target_company) then
            rank = pst_get_company_rank(i-1)
            if (rank .eq. 0) then
               master = i-1
               go to 50
            endif
         endif
      enddo
   50 continue

      call load_pre_defined_routines()

      if (master .eq. me) then

c---------------------------------------------------------------------------
c   Compile the SIAL code (or read an object file) for this company.
c---------------------------------------------------------------------------

        sial_file = company_sial_prog(i)
        if (sial_file(1:1) .eq. ' ') go to 100 ! no program

        n = str_trimlen(sial_file)
        if (sial_file(n-3:n) .eq. '.sio') then
           call read_tables(sial_file(1:n) // char(0))
        else if (sial_file(n-4:n) .eq. '.sial') then
            print *,'Error: Cannot use .sial files.'
            call abort_job()
c           ierr = interpretface(sial_file(1:n) // char(0) )
c           if (ierr .ne. 0) then
c              print *,'Error: SIAL compile errors.'
c              call abort_job()
c           endif
        else
           print *,'Error: Invalid file format for SIAL program.'
           print *,'File extension must be .sial or .sio'
           print *,'SIAL_FILE = ',sial_file(1:n)
           call abort_job()
        endif

         msg(1) = nindex_table_sip
         msg(2) = narray_table_sip
         msg(3) = noptable_sip
         msg(4) = nscalar_table_sip
         call mpi_bcast(msg, 4, MPI_INTEGER, 0,
     *               company_comm, ierr)
      else
         call mpi_bcast(msg, 4, MPI_INTEGER, 0,
     *               company_comm, ierr)

c------------------------------------------------------------------------------
c   Allocate memory in which to recv the tables.
c-----------------------------------------------------------------------------

         ierr = create_index_table(msg(1))
         if (ierr .ne. 0) then
            print *,'Task ',me,' Error: Cannot create index_table, ',
     *           'nindex_table = ',msg(1)
            call c_exit(1)
         endif

         ierr = create_array_table(msg(2))
         if (ierr .ne. 0) then
            print *,'Task ',me,' Error: Cannot create array_table, ',
     *           'narray_table = ',msg(2)
            call c_exit(1)
         endif

         ierr = create_optable(msg(3))
         if (ierr .ne. 0) then
            print *,'Task ',me,' Error: Cannot create optable, ',
     *           'narray_table = ',msg(3)
            call c_exit(1)
         endif

         nindex_table_sip = msg(1)
         narray_table_sip = msg(2)
         noptable_sip     = msg(3)
         nscalar_table_sip = msg(4)  
      endif      ! master .eq. me

c-----------------------------------------------------------------------------
c   Broadcast the 4 tables.
c   Index table:
c----------------------------------------------------------------------------

      n = nindex_table_sip*lindex_table_entry
      ind = get_index_from_base(index_table_base, ix, 1)
      call mpi_bcast(ix(ind), n, MPI_INTEGER, 0, company_comm, ierr)

c---------------------------------------------------------------------------
c   Array table
c---------------------------------------------------------------------------

      n = narray_table_sip*larray_table_entry
      ind = get_index_from_base(array_table_base, ix, 1)
      call mpi_bcast(ix(ind), n, MPI_INTEGER, 0, company_comm, ierr)
      
c---------------------------------------------------------------------------
c   Operation table
c---------------------------------------------------------------------------

      n = noptable_sip*loptable_entry
      ind = get_index_from_base(optable_base, ix, 1)
      call mpi_bcast(ix(ind), n, MPI_INTEGER, 0, company_comm, ierr)
      
c---------------------------------------------------------------------------
c   Scalar table
c---------------------------------------------------------------------------

      n = nscalar_table_sip
      ind = get_index_from_base(scalar_table_base, x, 2)
      call mpi_bcast(x(ind), n, MPI_DOUBLE_PRECISION, 0, 
     *                company_comm, ierr)
      

  100 continue   

      return
      end
