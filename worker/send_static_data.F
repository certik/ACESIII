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
      subroutine send_static_data(need_predef, npre_defined)
c---------------------------------------------------------------------------
c   Post the mpi_isends needed to distribute the static data to each worker.
c---------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'proto_events.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer npre_defined
      logical need_predef(npre_defined)

      double precision dbuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
#else
      common dbuf
#endif
      integer master, pst_get_master, pst_get_company
      integer i
      integer mpierr
      
#ifdef ALTIX
      dptr = dshptr
#endif
      master = pst_get_master()
      if (me .eq. master) then

c--------------------------------------------------------------------------
c   Check that the io_company_id is different from the company id.
c--------------------------------------------------------------------------

         do i = 1, nprocs
            scfa_req(i) = mpi_request_null
            scfb_req(i) = mpi_request_null
            epsa_req(i) = mpi_request_null
            epsb_req(i) = mpi_request_null
            focka_req(i) = mpi_request_null
            fockb_req(i) = mpi_request_null
            if (pst_get_company(i-1) .ne. io_company_id) then 
               if (need_predef(1) .or.
     *             need_predef(2)) 
     *            call mpi_isend(dbuf(iscfa), nscfa,
     *                           mpi_double_precision, i-1,
     *                           scfa_coeff_request_event,
     *                           mpi_comm_world, scfa_req(i), mpierr)

               if (need_predef(3))
     *            call  mpi_isend(dbuf(iscfb), nscfb,
     *                           mpi_double_precision, i-1,
     *                           scfb_coeff_request_event,
     *                           mpi_comm_world, scfb_req(i), mpierr)

               if (iuhf .eq. 2) then
                  if (need_predef(9))
     *               call mpi_isend(dbuf(ifockrohfa), nscfa,
     *                           mpi_double_precision, i-1,
     *                           focka_coeff_request_event,
     *                           mpi_comm_world, focka_req(i), mpierr)
              
                  if (need_predef(10))
     *               call  mpi_isend(dbuf(ifockrohfb), nscfb,
     *                           mpi_double_precision, i-1,
     *                           fockb_coeff_request_event,
     *                           mpi_comm_world, fockb_req(i), mpierr)
               endif

               call  mpi_isend(dbuf(iepsa), nepsa,
     *                           mpi_double_precision, i-1,
     *                           epsa_coeff_request_event,
     *                           mpi_comm_world, epsa_req(i), mpierr)
               if (nepsb .gt. 0)
     *            call  mpi_isend(dbuf(iepsb), nepsb,
     *                           mpi_double_precision, i-1,
     *                           epsb_coeff_request_event,
     *                           mpi_comm_world, epsb_req(i), mpierr)
            endif
         enddo
      endif

      return
      end
