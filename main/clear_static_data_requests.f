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
      subroutine clear_static_data_requests()
c----------------------------------------------------------------------------
c   Remove any outstanding mpi requests from sending the static data.
c----------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'proto_events.h'

      integer pst_get_company
      integer i, mpierr
      integer status(MPI_STATUS_SIZE)
      logical flag

      if (me .eq. 0) then
         do i = 1, nprocs
            if (pst_get_company(i-1) .ne. io_company_id) then 
               if (scfa_req(i) .ne. mpi_request_null) then
                  call mpi_test(scfa_req(i), flag, status, mpierr)
                  if (.not. flag) 
     *               call mpi_request_free(scfa_req(i), mpierr)
               endif

               if (scfb_req(i) .ne. mpi_request_null) then
                   call mpi_test(scfb_req(i), flag, status, mpierr)
                   if (.not. flag) 
     *                call mpi_request_free(scfb_req(i), mpierr)
               endif

               if (focka_req(i) .ne. mpi_request_null) then
                   call mpi_test(focka_req(i), flag, status, mpierr)
                   if (.not. flag)  
     *                call mpi_request_free(focka_req(i), mpierr)
               endif

               if (fockb_req(i) .ne. mpi_request_null) then
                   call mpi_test(fockb_req(i), flag, status, mpierr)
                   if (.not. flag)
     *                call mpi_request_free(fockb_req(i), mpierr)
               endif

               if (epsa_req(i) .ne. mpi_request_null) then
                   call mpi_test(epsa_req(i), flag, status, mpierr)
                   if (.not. flag)
     *                call mpi_request_free(epsa_req(i), mpierr)
               endif

               if (epsb_req(i) .ne. mpi_request_null) then
                  call mpi_test(epsb_req(i), flag, status, mpierr)
                   if (.not. flag)
     *                call mpi_request_free(epsb_req(i), mpierr)
               endif
            endif
         enddo
      endif

      return
      end
