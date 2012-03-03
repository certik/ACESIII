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
      subroutine sip_integral_init()
c--------------------------------------------------------------------
c
c   sip_integral_init: Handles initialization of worker tasks for
c                         the direct integral computation used by 
c                         the sip implementation.
c 
c---------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'proto_defines.h'
      include 'proto_events.h'
      include 'mem.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'dbugcom.h'
      include 'erd_dbug.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer me, ierr, my_company, company_comm, company_size
      integer iBytes, i0int
      integer i, j, ncomp
      integer master
      integer itemp, ztemp
      integer nalpha
      integer ihrr, zhrr 

      integer nw
      integer*8 nw8
      integer*8 c_loc64
      integer*8 get_index_from_base
      integer*8 addr
      integer status(MPI_STATUS_SIZE)

      integer pst_get_role
      integer pst_get_my_company
      integer pst_get_company
      integer pst_get_mem_scratch
      integer pst_get_master
      integer pst_get_company_comm

      logical configure
      logical spherical
      
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (iptr, ibuf)
#else
      equivalence (ibuf, dbuf)  
      common dbuf    ! blank common for scratch.
#endif

      logical do_erd_scratch_calc
      save do_erd_scratch_calc

      data do_erd_scratch_calc/.true./

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      my_company = pst_get_my_company()
      company_comm = pst_get_company_comm(me) 
      call mpi_comm_size(company_comm, company_size, ierr)
      master = pst_get_master()
         
#ifdef ALTIX
      dptr = dshptr     ! set dbuf pointer to shared mem.
      iptr = ishptr     ! set ibuf pointer to shared mem.
#endif

c------------------------------------------------------------------------
c   Initialization code for integral company.
c------------------------------------------------------------------------

      configure = .true.
      if (intpkg .eq. flocke_package) then
         if (master .eq. me .and. 
     *       .not. master_is_worker) then
            configure = .false.
         else
            nw = nalloc_pkg
            if (nw/2*2 .ne. nw) nw = nw + 1   
            call mem_alloc(dbuf, nw, bytes_per_double, icoord, .true.,
     *                     ierr)
            if (ierr .ne. 0) then
               print *,'Task ',me,' mem_alloc error: basis functions ',
     *                 nw
               call abort_job()
            endif

            nw = nbasis
            if (nw/2*2 .ne. nw) nw = nw + 1   
            call mem_alloc(dbuf, nw, bytes_per_double, iscale_fac, 
     *                     .true., ierr)
            if (ierr .ne. 0) then
               print *,'Task ',me,' mem_alloc error: scale factor',
     *                 nw
               call abort_job()
            endif

            call mem_alloc(ibuf, nw, intsize, ierdind, 
     *                     .true., ierr)
            if (ierr .ne. 0) then
               print *,'Task ',me,' mem_alloc error: ierdind',
     *                 nw
               call abort_job()
            endif

c-------------------------------------------------------------------------
c   Determine max. memory for "ccbeg", "ccend" arrays used in ERD package.
c-------------------------------------------------------------------------

            nw = 0
            do i = 1, nshells
               nw = nw + ncfps(i)
            enddo
            if (nw/2*2 .ne. nw) then
               nw = nw + 1
            else
               nw = nw + 2
            endif
 
            call mem_alloc(ibuf, nw, intsize, iccbeg,
     *                     .true., ierr)
            if (ierr .ne. 0) then
               print *,'Task ',me,' mem_alloc error: iccbeg',
     *                 nw
               call abort_job()
            endif

            call mem_alloc(ibuf, nw, intsize, iccend,
     *                     .true., ierr)
            if (ierr .ne. 0) then
               print *,'Task ',me,' mem_alloc error: iccend',
     *                 nw
               call abort_job()
            endif
         endif
      endif

      if (intpkg .eq. flocke_package) then

c-------------------------------------------------------------------------
c   Recv the obtain the basis functions from the master.
c-------------------------------------------------------------------------

         master = pst_get_master()
         if (master .ne. me) then
            call mpi_recv(dbuf(icoord), nalloc_pkg, 
     *                  mpi_double_precision, master,
     *                  basis_function_request_event, 
     *                  company_comm, status, ierr)  
         else
            icoord = master_icoord
         endif

         call mpi_barrier(company_comm, ierr)

         nalpha = 0
         do i = 1, nshells
            nalpha = nalpha + npfps(i)
         enddo

         ialpha  = icoord + 3*nshells
         ipcoeff = ialpha + nalpha

         if (configure) then
            call aces_to_erd(nfps, ivangmom, nshells, ispherical,
     *                       ibuf(ierdind), dbuf(iscale_fac)) 

            call setup_ccbeg(dbuf(ialpha), ixalpha, dbuf(ipcoeff), 
     *                       ixpcoef, ncfps, npfps, nshells, 
     *                       ibuf(iccbeg), ibuf(iccend), indx_cc)
         endif

         spherical = (ispherical .eq. 1) 

         if (do_erd_scratch_calc) then
            do_erd_scratch_calc = .false.   ! only 1 time EVER!

            if (fast_erd_memcalc) then
               call fast_erd_scratch_mem_calc(nshells, ivangmom, 
     *                          ncfps, npfps,
     *                          atom, dbuf(icoord), dbuf(ialpha),
     *                          dbuf(ipcoeff), ixalpha, ixpcoef,
     *                          ibuf(iccbeg), ibuf(iccend), indx_cc,
     *                          spherical, ncenters, intmax, zmax)
            else
               call erd_scratch_mem_calc(nshells, ivangmom, 
     *                          ncfps, npfps,
     *                          atom, dbuf(icoord), dbuf(ialpha),
     *                          dbuf(ipcoeff), ixalpha, ixpcoef,
     *                          ibuf(iccbeg), ibuf(iccend), indx_cc, 
     *                          spherical, ncenters, 
     *                          calc_2der_integrals,
     *                          intmax, zmax, me, company_size)
            endif

c-------------------------------------------------------------------------
c   Add correction for memory needed in the HRR transformation.
c-------------------------------------------------------------------------

            call erd__memory_hrr_correction(ivAngMom, nshells, 
     *                  spherical, ihrr, zhrr)
            intmax = intmax + ihrr
            zmax   = zmax + zhrr  
         endif

         itemp = intmax
         ztemp = zmax
         call mpi_barrier(company_comm, ierr)
         call mpi_reduce(itemp, intmax, 1,
     *          mpi_integer, mpi_max, 0,
     *          company_comm, ierr)
         call mpi_bcast(intmax, 1, MPI_INTEGER, 0,
     *               company_comm, ierr)

         call mpi_reduce(ztemp, zmax, 1,
     *                mpi_integer, mpi_max, 0,
     *                company_comm, ierr)
         call mpi_bcast(zmax, 1, MPI_INTEGER, 0,
     *               company_comm, ierr)

         if (me .eq. 0 .and. dbg) 
     *       print *,'Integral memory: intmax ',
     *        intmax, ' zmax ',zmax  

c---------------------------------------------------------------------------
c   Add the integral scratch buffers to the memory configuration, and
c   place it in context.
c---------------------------------------------------------------------------

         nw        = zmax
         call mem_alloc(dbuf, nw, bytes_per_double, dpkgscr,
     *                     .true., ierr)
         if (ierr .ne. 0) then
            print *,'Task ',me,' mem_alloc error: dpkgscr',
     *                 nw
            call abort_job()
         endif

         nw        = intmax
         call mem_alloc(ibuf, nw, intsize, ipkgscr,
     *                     .true., ierr)
         if (ierr .ne. 0) then
            print *,'Task ',me,' mem_alloc error: ipkgscr',
     *                 nw
            call abort_job()
         endif
      endif

         max_int_erd = intmax
         max_dp_erd  = zmax
         erd_base_int_addr = c_loc64(ibuf, ipkgscr, intsize)
         erd_base_dp_addr  = c_loc64(dbuf, dpkgscr, bytes_per_double)
         nw8 = zmax
         erd_dp_bounds = c_loc64(dbuf(dpkgscr), nw8, bytes_per_double)
         nw8 = intmax  
         erd_int_bounds = c_loc64(ibuf(ipkgscr), nw8, intsize)
         if (erd_dp_bounds .gt. erd_base_int_addr) then
            print *,'Task ',me,' INTEGRAL MEM ALLOC PROBLEM'
            print *,'Task ',me,' Base of dp mem: ',erd_base_dp_addr
            print *,'Task ',me,' End of dp mem = ',erd_dp_bounds
            print *,'Task ',me,' max_dp_erd    = ',max_dp_erd
            print *,'Task ',me,' Base of int mem = ',erd_base_int_addr
            print *,'Task ',me,' End of int mem = ',erd_int_bounds
            print *,'Task ',me,' max_int_erd    = ',max_int_erd
            call abort_job()
         endif  

c--------------------------------------------------------------------------
c   Initialize flag for computation of 1-electron integrals.
c--------------------------------------------------------------------------

      compute_1e_integrals = .true.

      return
      end

