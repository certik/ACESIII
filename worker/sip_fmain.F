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
      subroutine sip_fmain_init(ierr_return, ncompany_workers_min)
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'int_gen_parms.h'
      include 'symbolic_constants.h'
      include 'sip_tables.h'
      include 'machine_types.h'
      include 'timerz.h'
      include 'proto_events.h' 
      include 'parallel_info.h'
      include 'dropmo.h'
      include 'sial_config_params.h' 
      include 'dbugcom.h'
      include 'trace.h'
      include 'hess.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

#ifdef HP
      integer*8 lenshmem
#else
      integer lenshmem
#endif
      integer mx_ao_seg_size, mx_occ_seg_size, mx_virt_seg_size
      integer nocc_orbitals, nvirt_orbitals 
      integer naocc_orbitals, nbocc_orbitals
      integer ierr_return

      integer i, ii, j, n, flag, ierr, istat
      integer mpierr
      integer blocksize
      integer*8 get_index_from_base
      integer*8 c_loc64, ixx
      integer*8 i8temp, nxthess
      integer str_trimlen
      integer status(mpi_status_size)
      integer*8 i8size
      integer size
      integer igrad

      integer*8 iopt, iindex, iarray, iscalar
      integer*8 iaddr_table
      integer*8 isegbase, isegtable, ibmap_table
      integer*8 iexp_index_table
      integer*8 iblkmgr, malloced_previously
      integer*8 iblk_header, iblkmgr_data, ifree_stack, iblk_addr
      integer*8 iblk_in_use
      integer*8 addr, addr_comm_buffer
      integer noptable, nindex_table, narray_table, nscalar_table
      integer nsegment_table
      integer nblk_overhead, nblocks_total

      integer optable(1), array_table(1), index_table(1)
      integer exp_index_table(1)
      integer segment_table(1)
      integer*8 address_table(2)
      integer proctab(200)
      double precision scalar_table(1)
      integer block_map_table(1)
      integer nblkmgr, nblkmgr_ints
      double precision blkmgr_scr(1)
      integer iblkmgr_scr(1)
      integer niocompany
      integer range, maxseg, nwork, nwork2
      integer*8 nw_total
      integer*8 irange1, irange2
      integer itotenerg, comm_timer, master, nsial_lines
      integer iscfeneg, istabvalue
      character*40 timerdesc 
      character*10 array_labels   

      integer company_comm
      integer pst_get_company_comm, pst_get_company_rank
      integer pst_get_company, pst_get_master

      integer nstacks
      integer stack_blocksizes(max_stacks)
      integer stack_blocks(max_stacks)

      integer npre_defined
      parameter (npre_defined = 19)
      logical need_predef(npre_defined)

      double precision x(2)
      integer double_ratio
      integer nblock_map_table
      integer nw
      integer ncompany_workers, ncompany_workers_min
      integer*8 ixca, ixcb, ixfa, ixfb
#ifdef ALTIX
      pointer (dptr, x)
      pointer (ptr_optable, optable)
      pointer (ptr_index_table, index_table)
      pointer (ptr_array_table, array_table)
      pointer (ptr_scalar_table, scalar_table)
      pointer (ptr_block_map, block_map_table)
      pointer (ptr_segment_table, segment_table)
      pointer (ptr_address_table, address_table)
      pointer (ptr_blkmgr_scr, blkmgr_scr)
#else
      common x
#endif
      double precision mygrad(3*ncenters) 

      logical debug, validate
      logical atom_based
      save

      call mpi_comm_rank(mpi_comm_world, me, mpierr)
      call mpi_comm_size(mpi_comm_world, nprocs, mpierr)
      company_comm    = pst_get_company_comm(me)
      my_company_rank = pst_get_company_rank(me)
      call mpi_comm_size(company_comm, ncompany_workers, mpierr)

      ierr_return = 0

#ifdef ALTIX
      dptr = dshptr
      ptr_optable      = ishptr
      ptr_index_table  = ishptr
      ptr_array_table  = ishptr
      ptr_scalar_table = dshptr
      ptr_block_map    = ishptr
      ptr_segment_table=ishptr
      ptr_address_table =i8shptr
      ptr_blkmgr_scr   = dshptr
#endif

c--------------------------------------------------------------------------
c   Set up the server for this processor.
c   First, count the number of processes in the IOCOMPANY.
c   Second, determine the server for this process based on the rank of 
c      the process within its company, and the number of prcesses in
c      the IOCOMPANY.
c--------------------------------------------------------------------------

      niocompany = 0
      if (io_company_id .gt. 0) then
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) 
     *          niocompany = niocompany + 1
         enddo
      endif

      iopt   = get_index_from_base(optable_base, optable, 1)
      iindex = get_index_from_base(index_table_base,
     *     index_table, 1)
      iarray = get_index_from_base(array_table_base, array_table, 1)

      iscalar = get_index_from_base(scalar_table_base, scalar_table, 2)

      noptable = noptable_sip
      nindex_table = nindex_table_sip
      narray_table = narray_table_sip
      nscalar_table = nscalar_table_sip

c---------------------------------------------------------------------
c   Set up the segment table.
c---------------------------------------------------------------------

      if (vvvi_stack) then
         nstacks = 9
      else
         nstacks = 7
      endif

      nocc_orbitals = noccupied
      naocc_orbitals = nalpha_occupied
      if (nbeta_occupied .ne. 0) then
         nbocc_orbitals = nbeta_occupied
      else
         nbocc_orbitals = nocc_orbitals
      endif

       mx_ao_seg_size = sip_mx_ao_segsize
       mx_occ_seg_size = sip_mx_occ_segsize
       mx_virt_seg_size = sip_mx_virt_segsize

c--------------------------------------------------------------------------
c   If dropmo is used, create a copy of the index_table, for later use
c   in expanding the basis.
c--------------------------------------------------------------------------

       if (.not. ignore_dropmo) then
          call mem_alloc(exp_index_table, 
     *            nindex_table*lindex_table_entry, intsize,
     *            iexp_index_table, .true., ierr)

          eindex_table_base = c_loc64(exp_index_table, 
     *                             iexp_index_table, intsize)

          do i = 1, nindex_table*lindex_table_entry
             exp_index_table(iexp_index_table+i-1) = 
     *             index_table(iindex+i-1)
          enddo   
       endif

c---------------------------------------------------------------------------
c   Determine whether we use atom_based AO segmentation or the normal
c   algorithm.
c---------------------------------------------------------------------------

       atom_based = .false.
       if (aat_database .ne. ' ') atom_based = .true.

c---------------------------------------------------------------------------
c   Build the segment table.
c---------------------------------------------------------------------------

       call build_segment_table(index_table(iindex), nindex_table,
     *                        end_nfps, nshells, mx_ao_seg_size,
     *                        mx_occ_seg_size, mx_virt_seg_size,
     *                        sip_sub_segsize, sip_sub_occ_segsize,
     *                        sip_sub_virt_segsize, sip_sub_ao_segsize,
     *                        nocc_orbitals, naocc_orbitals,
     *                        nbocc_orbitals, 
     *                        nvirt_orbitals, nalpha_virtual,
     *                        nbeta_virtual, 
     *                        isegbase,nsegment_table,
     *                        stack_blocksizes, nstacks,
     *                        array_table(iarray), narray_table,
     *                        atom_based, atom)

       nsegment_table_sip = nsegment_table  ! store in common

c----------------------------------------------------------------------------
c   Add additional symbolic constants from parameter data.
c----------------------------------------------------------------------------

       symbolic_constant_table(38) = cc_iter
       symbolic_constant_table(39) = cc_hist
       symbolic_constant_table(40) = cc_beg
       symbolic_constant_table(41) = scf_iter
       symbolic_constant_table(42) = scf_hist
       symbolic_constant_table(43) = scf_beg
       symbolic_constant_table(44) = ncenters
       symbolic_constant_table(45) = itrips
       symbolic_constant_table(46) = itripe
       symbolic_constant_table(47) = ihess1
       symbolic_constant_table(48) = ihess2
       symbolic_constant_table(49) = jhess1
       symbolic_constant_table(50) = jhess2
       symbolic_constant_table(51) = subb
       symbolic_constant_table(52) = sube
       symbolic_constant_table(53) = sip_sub_segsize
       symbolic_constant_table(54) = sip_sub_occ_segsize
       symbolic_constant_table(55) = sip_sub_virt_segsize
       symbolic_constant_table(56) = sip_sub_ao_segsize

       isegtable = get_index_from_base(isegbase, segment_table, 1)
       segment_table_base = isegbase
            
c------------------------------------------------------------------------
c   Allocate enough memory for the address_table.
c------------------------------------------------------------------------

      size = loc(address_table(2)) - loc(address_table(1))

      call mem_alloc(address_table, narray_table, size, iaddr_table,
     *               .true., ierr)
      ixx = 1
      address_table_base = c_loc64(address_table(iaddr_table), ixx, 8)
      proctab_base       = c_loc64(proctab, ixx, intsize)

c--------------------------------------------------------------------------
c   Determine the amount of memory needed for the block_map_table.
c--------------------------------------------------------------------------

      call block_map_table_setup(array_table(iarray), narray_table,
     *                     index_table(iindex), nindex_table,
     *                     segment_table(isegbase), nsegment_table,
     *                     nblock_map_table)

c----------------------------------------------------------------------------
c   Allocate the memory for the block_map_table.
c----------------------------------------------------------------------------

      nw = nblock_map_table * lblock_map_entry
      call mem_alloc(block_map_table, nw, intsize, 
     *               ibmap_table,
     *               .true., ierr)

      ixx = 1
      bmap_table_base = c_loc64(block_map_table, ibmap_table, 
     *                       intsize)

c------------------------------------------------------------------------
c   Initialize all the tables.
c------------------------------------------------------------------------

       call vtdemo_init(optable(iopt), noptable, 
     *                    array_table(iarray), narray_table,
     *                    index_table(iindex), nindex_table, 
     *                    segment_table(isegtable), nsegment_table,
     *                    scalar_table(iscalar), nscalar_table, 
     *                    block_map_table(ibmap_table),
     *                    nblock_map_table, proctab, 
     *                    address_table(iaddr_table), blocksize,
     *                    end_nfps, nshells, scf_energy, totenerg,
     *                    damp_init, cc_conv, scf_conv, stabvalue,
     *
     *                    excite, eom_tol, eom_roots, ! Watson Added
     *
     *                    io_company_id, niocompany, need_predef,
     *                    npre_defined, dryrun)
 
c------------------------------------------------------------------------
c   Determine block size for the thread buffer(s).
c------------------------------------------------------------------------

      blocksize = 0
      do i = 1, nstacks
         blocksize = max(stack_blocksizes(i), blocksize)
      enddo
      
c--------------------------------------------------------------------------
c   Allocate the remaining memory for use as data blocks.
c--------------------------------------------------------------------------

      i8temp = maxmem * 1024
      i8size = i8temp * 1024
      call mem_alloc_query(malloced_previously)
      i8size    = i8size - malloced_previously
      i8size    = i8size / bytes_per_double
      nblkmgr = (0.95*i8size) 
      call mem_alloc(blkmgr_scr, nblkmgr, bytes_per_double, 
     *                  iblkmgr, .true., ierr)
      if (ierr .ne. 0) then
         print *,'Task ',me,' mem_alloc error in sip_fmain_init ',
     *        nblkmgr
         call abort_job()
      endif

      i8size = nblkmgr
      i8size = i8size * bytes_per_double
      if (me .eq. 0) print *,'Allocated ',i8size,
     *        ' bytes for blkmgr.'

c---------------------------------------------------------------------------
c   Set up the addresses of communication buffers.
c---------------------------------------------------------------------------

      addr_comm_buffer = c_loc64(blkmgr_scr,iblkmgr,
     *                                    bytes_per_double)

      ithread = get_index_from_base(addr_comm_buffer, x, 2)
      if (nprocs .gt. 1) then
         nwork   = nworkthread * blocksize 
      else
         nwork = 0    ! thread server will not be active.
      endif

c-------------------------------------------------------------------------
c   Do a dryrun to determine the proper distribution of the stack blocks.
c-------------------------------------------------------------------------

       nblk_overhead = lblk_header + 2 + 8/intsize  ! overhead per block

       call stack_distribution(optable(iopt), noptable, 
     *                    array_table(iarray), narray_table,
     *                    index_table(iindex), nindex_table, 
     *                    segment_table(isegtable), nsegment_table,
     *                    scalar_table(iscalar), nscalar_table, 
     *                    block_map_table(ibmap_table),
     *                    nblock_map_table, proctab, 
     *                    blkmgr_scr(iblkmgr+nwork), nblkmgr-nwork,
     *                    nblk_overhead,
     *                    stack_blocksizes, nstacks, stack_blocks,
     *                    stack_algorithm_type, ncompany_workers,
     *                    niocompany, ierr, .true.)
      if (ierr .ne. 0) then
         
c--------------------------------------------------------------------------
c   Attempt to determine the number of processors on which this job 
c   will fit.
c--------------------------------------------------------------------------

  100    continue
         ncompany_workers = ncompany_workers + 1
         call refill_block_map_table(ncompany_workers, 
     *                    block_map_table(ibmap_table),
     *                    nblock_map_table, 
     *                    array_table(iarray), narray_table,
     *                    index_table(iindex), nindex_table)
         call stack_distribution(optable(iopt), noptable,
     *                    array_table(iarray), narray_table,
     *                    index_table(iindex), nindex_table,
     *                    segment_table(isegtable), nsegment_table,
     *                    scalar_table(iscalar), nscalar_table,
     *                    block_map_table(ibmap_table),
     *                    nblock_map_table, proctab,
     *                    blkmgr_scr(iblkmgr+nwork), nblkmgr-nwork,
     *                    nblk_overhead,
     *                    stack_blocksizes, nstacks, stack_blocks,
     *                    stack_algorithm_type, ncompany_workers,
     *                    niocompany, istat, .false.)
         if (istat .ne. 0 .and. 
     *       ncompany_workers .lt. max_procsx) go to 100
      endif

      call mpi_barrier(company_comm, mpierr)
      call mpi_reduce(ncompany_workers, ncompany_workers_min, 1,
     *                mpi_integer, mpi_max, 0,
     *                company_comm, mpierr)

      ierr_return = ierr

      if (my_company_rank .eq. 0 .and. ierr .ne. 0) then
         print *,'DRYRUN ERROR: Step will abort, ',
     *       'cannot fit in memory'
         print *,'This step requires a minimum of ',
     *      ncompany_workers_min,
     *      ' workers along with ',niocompany,' servers to run.'
         return
      endif

c-------------------------------------------------------------------------
c   Create the server table and send a copy to each server.
c-------------------------------------------------------------------------

      if (my_company_rank .eq. 0 .and. .not. dryrun) then
         i8size = nblkmgr
         i8size = i8size * bytes_per_double
         nblkmgr_ints = i8size / intsize

         i8size = niocompany * nstacks
         i8size = i8size * intsize
         i8size = (i8size + bytes_per_double -1)/ bytes_per_double
         nwork2 = i8size 

         nblkmgr_ints = nblkmgr_ints - 
     *                  (nwork2*bytes_per_double) / intsize
         
         call distribute_server_info(array_table(iarray), narray_table,
     *                     index_table(iindex), nindex_table,
     *                     segment_table(isegtable), nsegment_table,
     *                     block_map_table(ibmap_table),
     *                     nblock_map_table, io_company_id,
     *                     niocompany, blkmgr_scr(iblkmgr+nwork2), 
     *                     nblkmgr_ints, stack_blocksizes, nstacks,
     *                     blkmgr_scr(iblkmgr))
      endif

c------------------------------------------------------------------------
c   Set up the following dynamically allocated arrays:
c      blk_header - For data required by blkmgr to manage each block.
c      free_stack - Array containing the free blocks for each stack.
c      blk_addr   - Array continuing the address of each data block.
c   These arrays will be initialized by blkmgr_init, but we must 
c   allocate them here.
c------------------------------------------------------------------------

      nblocks_total = 0
      nw_total      = 0
      do i = 1, nstacks
         nblocks_total = nblocks_total + stack_blocks(i)
         nw_total = nw_total + stack_blocks(i)*(stack_blocksizes(i) +
     *                2) 
      enddo

      iblkmgr_data = iblkmgr + nwork
      nblkmgr      = nblkmgr - nwork

      nblkmgr = nblkmgr - (nblk_overhead * nblocks_total *
     *                      intsize) / bytes_per_double

      ixx = iblkmgr_data + nblkmgr   ! point to end of block area.
      addr = c_loc64(blkmgr_scr, ixx, bytes_per_double)
      iblk_header = get_index_from_base(addr, iblkmgr_scr, 1)
      ifree_stack = iblk_header + lblk_header * nblocks_total
      iblk_in_use = ifree_stack + nblocks_total
      iblk_addr   = iblk_in_use + nblocks_total 

c-------------------------------------------------------------------------
c   Force the blk_addr table to come out on a 8-byte boundary.
c-------------------------------------------------------------------------

      addr = c_loc64(iblkmgr_scr, iblk_addr, intsize)
      if (addr/8 * 8 .ne. addr) then
         addr = (addr + 7)/8 * 8
         iblk_addr = get_index_from_base(addr, iblkmgr_scr, 1)
      endif
 
c------------------------------------------------------------------------
c   Initialize the block manager.
c------------------------------------------------------------------------

       call blkmgr_init(blkmgr_scr(iblkmgr_data), nblkmgr, 
     *                  nblocks_total,
     *                  iblkmgr_scr(iblk_header), 
     *                  iblkmgr_scr(ifree_stack),
     *                  iblkmgr_scr(iblk_addr),
     *                  iblkmgr_scr(iblk_in_use),
     *                  stack_blocksizes, stack_blocks, 
     *                  nstacks, company_comm)

      if (me .eq. 0) then
         call mem_alloc_query(malloced_previously)
         i8size = malloced_previously
         i8size = i8size / 1024
         i8size = i8size / 1024
         print *,'Total memory usage ',i8size,' MBytes.'
         print *,'Max. possible usage ',maxmem,' MBytes' 
         print *,'Total blocks used = ',nblocks_total
      endif 
      return

      entry sip_fmain()
      call mpi_comm_rank(mpi_comm_world, me, ierr)
      company_comm    = pst_get_company_comm(me)
      my_company_rank = pst_get_company_rank(me)

#ifdef ALTIX
      dptr = dshptr
      ptr_optable      = ishptr
      ptr_index_table  = ishptr
      ptr_array_table  = ishptr
      ptr_scalar_table = dshptr
      ptr_block_map    = ishptr
      ptr_segment_table=ishptr 
      ptr_address_table = i8shptr
      ptr_blkmgr_scr   = dshptr
#endif

c----------------------------------------------------------------------------
c   Set up table addresses.
c----------------------------------------------------------------------------

      ixx = 1
      iopt   = get_index_from_base(optable_base, optable, 1)
      iindex = get_index_from_base(index_table_base,
     *     index_table, 1)
      iarray = get_index_from_base(array_table_base, array_table, 1)

      iscalar = get_index_from_base(scalar_table_base, scalar_table, 2)
      iaddr_table = get_index_from_base(address_table_base, 
     *                                  address_table, 3)
      ibmap_table = get_index_from_base(bmap_table_base,
     *                                 block_map_table, 1)
      isegtable = get_index_from_base(segment_table_base, segment_table,
     *                                1)

      noptable = noptable_sip
      nindex_table = nindex_table_sip
      narray_table = narray_table_sip
      nscalar_table = nscalar_table_sip

      timerdesc = 'Block wait time'
      call register_timer(timerdesc, elapsed_time_timer, comm_timer)

      debug = .false.
      validate = .false.
      
      call set_sumz_tables(array_table(iarray), 
     *                     block_map_table(ibmap_table))

      call optable_loop(optable(iopt), noptable, array_table(iarray),
     *                narray_table,
     *                array_labels, index_table(iindex), nindex_table,
     *                segment_table(isegtable), nsegment_table,
     *                block_map_table(ibmap_table), nblock_map_table,
     *                scalar_table(iscalar), nscalar_table,
     *                proctab, address_table(iaddr_table), 
     *                debug, validate,
     *                company_comm, comm_timer)

      if (my_company_rank .eq. 0) then
         print *, '---------- Finished! ---------------'

         if (dbg) then
            print *,'Final state of scalar table:'
            do i = 1, nscalar_table
              print *,'scalar entry ',i,': ',
     *              scalar_table(iscalar+i-1)
            enddo
         endif

c-------------------------------------------------------------------------
c   Save TOTENERG value for later update of JOBARC file.
c-------------------------------------------------------------------------

         itotenerg = 
     *          array_table(iarray+4*larray_table_entry+
     *                       c_scalar_index-1)
         print *,'Total energy: ',
     *          scalar_table(iscalar+itotenerg-1)
         totenerg = scalar_table(iscalar+itotenerg-1)

c-------------------------------------------------------------------------
c   Save SCFENEG value for later update of JOBARC file.
c-------------------------------------------------------------------------

         iscfeneg =
     *          array_table(iarray+3*larray_table_entry+
     *                      c_scalar_index-1)
         print *,'SCF energy: ',
     *          scalar_table(iscalar+iscfeneg-1)
         scf_energy = scalar_table(iscalar+iscfeneg-1)

c-------------------------------------------------------------------------
c   Save STABVALUE value for later update of JOBARC file.
c-------------------------------------------------------------------------

         istabvalue =
     *          array_table(iarray+18*larray_table_entry+
     *                      c_scalar_index-1)
         print *,'HF Instability Value: ',
     *          scalar_table(iscalar+istabvalue-1)
         stabvalue = scalar_table(iscalar+istabvalue-1)
      endif

c-------------------------------------------------------------------------
c   Sum up the distributed contributions to the gradient.
c-------------------------------------------------------------------------

      call mpi_barrier(company_comm, ierr)
      call mpi_comm_size(company_comm, ncompany_workers, ierr)
      if (ncompany_workers .gt. 1) then
         do i = 1, 3*ncenters
            mygrad(i) = gradient_data(i)
         enddo

         call mpi_reduce(mygrad, gradient_data, 3*ncenters,
     *                mpi_double_precision, mpi_sum, 0, 
     *                company_comm, ierr)

c----------------------------------------------------------------------------
c   Final summation of the Hessian is performed in batches to save memory.
c   The mygrad array is used as a temporary buffer so the send and receive 
c   buffers of the mpi_reduce call do not overlap.
c----------------------------------------------------------------------------

         nxthess = ix_hess
         do j = 1, 3*ncenters
            do i = 1, 3*ncenters
               mygrad(i) = x(nxthess+i-1)
            enddo

            call mpi_reduce(mygrad, x(nxthess), 3*ncenters,
     *                mpi_double_precision, mpi_sum, 0,
     *                company_comm, ierr)
            nxthess = nxthess + 3*ncenters
         enddo
      endif

      if (my_company_rank .eq. 0) then
         print *,'Gradient data:'
         do i = 1, ncenters
            igrad = (i-1)*3 + 1
            print *,'Atom #',i,': ',gradient_data(igrad),
     *               gradient_data(igrad+1), gradient_data(igrad+2)
         enddo
     
c---------------------------------------------------------------------------
c   Get addresses of CA, CB, FOCK_A, and FOCK_B arrays. 
c   ONLY if the scf has just been performed. 
c---------------------------------------------------------------------------

         if (if_scf .or. init_scf) then
            ixca = 1
            ixcb = 1
            ixfa = 1
            ixfb = 1
            if (address_table(iaddr_table+1) .ne. 0) 
     *         ixca = get_index_from_base(address_table(iaddr_table+1),
     *                                     x,2) 
            if (address_table(iaddr_table+2) .ne. 0)
     *         ixcb = get_index_from_base(address_table(iaddr_table+2),
     *                                    x,2) 
            if (address_table(iaddr_table+8) .ne. 0)
     *         ixfa = get_index_from_base(address_table(iaddr_table+8),
     *                                    x,2)
            if (address_table(iaddr_table+9) .ne. 0) 
     *         ixfb = get_index_from_base(address_table(iaddr_table+9),
     *                                    x,2) 

             call write_eig_data(x(ixca), x(ixcb), x(ixfa), x(ixfb), 
     *                           nbasis, need_predef)
         endif 

      endif

      return
      end
