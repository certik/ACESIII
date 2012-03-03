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
      subroutine vtdemo_init(optable, noptable, 
     *                       array_table, narray_table,
     *                       index_table, nindex_table,
     *                       segment_table, nsegment_table,
     *                       scalar_table, nscalar_table,
     *                       block_map_table, nblock_map_table,
     *                       proctab, address_table,
     *                       blocksize, end_nfps, nshells, scf_energy,
     *                       totenerg, damp_init, cc_conv, scf_conv,
     *                       stabvalue,
     *
     *                       excite, eom_tol, eom_roots, ! Watson Added
     *
     *                       io_company_id, niocompany, need_predef,
     *                       npre_defined, dryrun)
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'proto_events.h'
      include 'timerz.h'
      include 'machine_types.h'
      include 'where_table.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer mx_proctab
      parameter (mx_proctab = 1000)

      integer*8 ixx
      integer nindex_table, narray_table, noptable, nsegment_table
      integer nscalar_table
      integer segment_table(lsegment_table_entry,nsegment_table)
      character*10 array_labels(narray_table)
      integer array_table(larray_table_entry, narray_table)
      integer*8 address_table(narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer optable(loptable_entry, noptable)
      integer nblock_map_table
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer proctab(2,mx_proctab)
      double precision scalar_table(nscalar_table)
      
      integer nshells
      integer end_nfps(nshells)
      integer segs(nshells)
      integer nsegs, first
      integer blocksize, nbasis
      integer next, next_server, next_worker
      integer proc_timer, target, key, key1, key2

      integer i, j, k, n
      integer nproctab
      integer me, nprocs, ierr
      integer seg
      integer nprocs_job
      integer size
      integer lineno, nsial_lines
      integer*8 ixmem

      integer*8 ixc, ixfock_a, ixfock_b
      integer*8 ixfockrohf_a, ixfockrohf_b
      integer*8 ixoed_nai, ixoed_kin, ixoed_ovl
      integer*8 nbytes_static
      integer iproc
      integer bseg_arr(mx_array_index)
      integer eseg_arr(mx_array_index)
      integer nseg_arr(mx_array_index)
      integer ind_arr(mx_array_index)
      integer type
      integer nind, nblks
      integer segknt
      integer get_subroutine_handle
      integer op(loptable_entry)
      integer opcode
      integer bseg, eseg, bval1, bval2, eval1, eval2
      integer itotenerg
      integer master, pst_get_master, predef_req
      integer iscfeneg
      integer iscalar
      integer istabvalue
      integer iexcite, ieom_tol, ieom_roots ! Watson Added
      integer timertype
      integer index, val1, val2
      integer occ_index, index_type, ind, nindex
      integer result
      integer iop, max_batch
      integer ind_pardo(mx_array_index)
      double precision scf_energy, totenerg, stabvalue
      double precision damp_init, cc_conv, scf_conv
      double precision excite, eom_tol, eom_roots ! Watson Added

      double precision xxx(2)
#ifdef ALTIX
      pointer (dptr, xxx)
#endif

      character*12 get_operation_desc 
      integer      get_operation_timer_type
      character*12 opdesc 
      character*40 timerdesc
      character*1 typechar

      integer my_company_rank, ncworkers
      integer pst_get_company
      integer pst_get_company_rank
      integer pst_get_company_comm
      integer company_comm, iproc_company_rank
      integer io_company_id, niocompany
      integer ioproc(niocompany)

      integer iseg(mx_array_index)
      integer lookup, block_map_lookup
      integer stat_key
      logical debug, validate, doit
      logical dryrun

      integer npre_defined
      integer setline(npre_defined)
      logical need_predef(npre_defined)
      logical i_used_it
      integer*8 c_loc64

      integer push_do, pop_do
      save opdesc

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      company_comm    = pst_get_company_comm(me)
      my_company_rank = pst_get_company_rank(me)
      call mpi_comm_size(company_comm, ncworkers, ierr)

#ifdef ALTIX
      dptr = dshptr
#endif

      nprocs_job = nprocs
      debug    = .false.
      validate = .false.

      size = bytes_per_double

c--------------------------------------------------------------------------
c   If an IOCOMPANY is used for served arrays, set up an array of process
c   ranks (relative to the MPI_COMM_WORLD communicator.
c--------------------------------------------------------------------------

      if (io_company_id .ne. 0) then
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) then
               iproc_company_rank = pst_get_company_rank(i-1)
               ioproc(iproc_company_rank+1) = i-1
            endif
         enddo
      endif

c---------------------------------------------------------------------------
c   Clear the address table.
c---------------------------------------------------------------------------

      do i = 1, narray_table
         address_table(i) = 0
         array_table(c_block_map,i) = 0
      enddo

c----------------------------------------------------------------------------
c   Find all pre-defined arrays that are used in specific instructions.
c----------------------------------------------------------------------------

      do i = 1, npre_defined
         need_predef(i) = .false.
         setline(i)     = 0
      enddo

      do i = 1, noptable
         if (optable(c_opcode,i) .lt. call_op .or.
     *       optable(c_opcode,i) .gt. sp_ldindex_op) then

            if (optable(c_op1_array,i) .gt. 0 .and.
     *          optable(c_op1_array,i) .le. npre_defined) then
               need_predef(optable(c_op1_array,i)) = .true.
               if (setline(optable(c_op1_array,i)) .eq. 0) 
     *           setline(optable(c_op1_array,i)) = optable(c_lineno,i)
            endif

            if (optable(c_op2_array,i) .gt. 0 .and.
     *          optable(c_op2_array,i) .le. npre_defined) then
               need_predef(optable(c_op2_array,i)) = .true.
               if (setline(optable(c_op2_array,i)) .eq. 0)
     *            setline(optable(c_op2_array,i)) = optable(c_lineno,i)
            endif

            if (optable(c_result_array,i) .gt. 0 .and.
     *          optable(c_result_array,i) .le. npre_defined) then
               need_predef(optable(c_result_array,i)) = .true.
               if (setline(optable(c_result_array,i)) .eq. 0)
     *            setline(optable(c_result_array,i)) = 
     *                              optable(c_lineno,i)
            endif
        endif 
      enddo

   10 continue
      if (me .eq. 0) then
         do i = 1, npre_defined
            if (need_predef(i)) 
     *         print *,'Static pre-defined array #',i,
     *             ' is first used on line ',setline(i)
         enddo

c----------------------------------------------------------------------------
c   Post the send for data needed to initialize static arrays from the
c   master to the other workers.
c----------------------------------------------------------------------------

c         if (.not. dryrun) 
c     *      call send_static_data(need_predef, npre_defined) 
      endif

c---------------------------------------------------------------------------
c   Set up a local array for C.
c---------------------------------------------------------------------------

      nbytes_static = 0
      if (need_predef(1) .or. 
     *    need_predef(2)) then   
         n = 1
         do j = 1, array_table(c_nindex,1)
            index = array_table(c_index_array1+j-1,1)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,1) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,1) = val2
            k = array_table(c_index_range2+j-1,1) - 
     *          array_table(c_index_range1+j-1,1) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixc, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating C array'
            call abort_job()
         endif

         address_table(1) = c_loc64(xxx, ixc, bytes_per_double)
         nbytes_static = nbytes_static + n * size

c---------------------------------------------------------------------------
c   Set up a local array for CA.
c---------------------------------------------------------------------------

         n = 1
         do j = 1, array_table(c_nindex,2)
            index = array_table(c_index_array1+j-1,2)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,2) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,2) = val2
            k = array_table(c_index_range2+j-1,2) -
     *          array_table(c_index_range1+j-1,2) + 1
            n = n * k
         enddo

         address_table(2) = address_table(1)

         op(c_result_array) = 1    ! C is array #1. 
         if (.not. dryrun) then 
           call init_c(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op)
         endif
      endif   ! need_predef

c---------------------------------------------------------------------------
c   Set up a local array for CB.
c---------------------------------------------------------------------------

      if (need_predef(3)) then
         n = 1
         do j = 1, array_table(c_nindex,3)
            index = array_table(c_index_array1+j-1,3)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,3) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,3) = val2
            k = array_table(c_index_range2+j-1,3) -
     *          array_table(c_index_range1+j-1,3) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixc, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating CB array'
            call abort_job()
         endif
         address_table(3) = c_loc64(xxx, ixc, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 3   ! CB is array #3
         if (.not. dryrun) then
            call init_c(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op)
         endif
      endif

      if (.not. dryrun) then
         call init_eps(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            op)
      endif

      iscfeneg = array_table(c_scalar_index, 4)   
      scalar_table(iscfeneg) = scf_energy
      itotenerg = array_table(c_scalar_index, 5)
      scalar_table(itotenerg) = totenerg

c---------------------------------------------------------------------------
c   Set up a local array for FOCK_A.  (Pre-defined array table entry 9)
c---------------------------------------------------------------------------

      if (need_predef(9)) then
         n = 1
         do j = 1, array_table(c_nindex,9)
            index = array_table(c_index_array1+j-1,9)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,9) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,9) = val2
            k = array_table(c_index_range2+j-1,9) -
     *          array_table(c_index_range1+j-1,9) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixfock_a, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating FOCKA array'
            call abort_job()
         endif
         address_table(9) = c_loc64(xxx, ixfock_a, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 9    ! FOCK_A is array #9.
         if (.not. dryrun) then 
            call init_fock(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op)
         endif
      endif

c---------------------------------------------------------------------------
c   Set up a local array for FOCK_B.  (Pre-defined array table entry 10)
c---------------------------------------------------------------------------

      if (need_predef(10)) then
         n = 1
         do j = 1, array_table(c_nindex,10)
            index = array_table(c_index_array1+j-1,10)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,10) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,10) = val2
            k = array_table(c_index_range2+j-1,10) -
     *          array_table(c_index_range1+j-1,10) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixfock_b, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating FOCKB array'
            call abort_job()
         endif
         address_table(10) = c_loc64(xxx, ixfock_b, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 10    ! FOCK_B is array #10.
         if (.not. dryrun) then
            call init_fock(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op)
         endif
      endif

c---------------------------------------------------------------------------
c   Set up a local array for oed_nai.  (Pre-defined array table entry 11)
c---------------------------------------------------------------------------

      if (need_predef(11)) then
         n = 1
         do j = 1, array_table(c_nindex,11)
            index = array_table(c_index_array1+j-1,11)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,11) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,11) = val2
            k = array_table(c_index_range2+j-1,11) -
     *          array_table(c_index_range1+j-1,11) + 1
            n = n * k
         enddo
   
         call mem_alloc(xxx, n, size, ixoed_nai, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating NAI integrals array'
            call abort_job()
         endif
         address_table(11) = c_loc64(xxx, ixoed_nai, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 11    ! ONE_NAI is array #11.
         if (.not. dryrun) then 
            call init_oed_integrals(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op, 1)
         endif
      endif

c---------------------------------------------------------------------------
c   Set up a local array for oed_kin.  (Pre-defined array table entry 12)
c---------------------------------------------------------------------------

      if (need_predef(12)) then
         n = 1
         do j = 1, array_table(c_nindex,12)
            index = array_table(c_index_array1+j-1,12)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,12) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,12) = val2
            k = array_table(c_index_range2+j-1,12) -
     *          array_table(c_index_range1+j-1,12) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixoed_kin, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating KINETIC integrals array'
            call abort_job()
         endif
         address_table(12) = c_loc64(xxx, ixoed_kin, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 12    ! ONE_KIN is array #12.
         if (.not. dryrun) then
            call init_oed_integrals(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op, 2)
         endif
      endif

c---------------------------------------------------------------------------
c   Set up a local array for oed_ovl.  (Pre-defined array table entry 13)
c---------------------------------------------------------------------------

      if (need_predef(13)) then
         n = 1
         do j = 1, array_table(c_nindex,13)
            index = array_table(c_index_array1+j-1,13)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,13) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,13) = val2
            k = array_table(c_index_range2+j-1,13) -
     *          array_table(c_index_range1+j-1,13) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixoed_ovl, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating OVERLAP integrals array'
            call abort_job()
         endif
         address_table(13) = c_loc64(xxx, ixoed_ovl, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 13    ! ONE_OVL is array #13.
         if (.not. dryrun) then
            call init_oed_integrals(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op, 3)
         endif
      endif

c---------------------------------------------------------------------------
c   Initialization of "damp" value.
c---------------------------------------------------------------------------

      iscalar = array_table(c_scalar_index, 14)
      scalar_table(iscalar) = damp_init

c---------------------------------------------------------------------------
c   Initialization of "cc_conv" value.
c---------------------------------------------------------------------------

      iscalar = array_table(c_scalar_index, 15)
      scalar_table(iscalar) = cc_conv

c---------------------------------------------------------------------------
c   Initialization of "scf_conv" value.
c---------------------------------------------------------------------------

      iscalar = array_table(c_scalar_index, 16)
      scalar_table(iscalar) = scf_conv

c---------------------------------------------------------------------------
c   Set up a local array for FOCKROHF_A.  (Pre-defined array table entry 17)
c---------------------------------------------------------------------------

      if (need_predef(17)) then
         n = 1
         do j = 1, array_table(c_nindex,17)
            index = array_table(c_index_array1+j-1,17)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,17) = val1
            call get_index_segment(index, eseg, segment_table,
     *         nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,17) = val2
            k = array_table(c_index_range2+j-1,17) -
     *          array_table(c_index_range1+j-1,17) + 1
            n = n * k
         enddo
         call mem_alloc(xxx, n, size, ixfockrohf_a, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating ROHFA array'
            call abort_job()
         endif
         address_table(17) = c_loc64(xxx, ixfockrohf_a, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 17    ! FOCKROHF_A is array #17.
         if (.not. dryrun) then
            call init_fockrohf(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op)
         endif    
      endif

c---------------------------------------------------------------------------
c   Set up a local array for FOCKROHF_B.  (Pre-defined array table entry 18)
c---------------------------------------------------------------------------

      if (need_predef(18)) then
         n = 1
         do j = 1, array_table(c_nindex,18)
            index = array_table(c_index_array1+j-1,18)
            bseg = index_table(c_bseg, index)
            eseg = index_table(c_eseg, index)
            call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range1+j-1,18) = val1
            call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            array_table(c_index_range2+j-1,18) = val2
            k = array_table(c_index_range2+j-1,18) -
     *          array_table(c_index_range1+j-1,18) + 1
            n = n * k
         enddo

         call mem_alloc(xxx, n, size, ixfockrohf_b, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Mem failure allocating ROHFB array'
            call abort_job()
         endif
         address_table(18) = c_loc64(xxx, ixfockrohf_b, size)
         nbytes_static = nbytes_static + n * size

         op(c_result_array) = 18    ! FOCKROHF_B is array #18.
         if (.not. dryrun) then
            call init_fockrohf(array_table, narray_table,
     *            index_table,
     *            nindex_table, segment_table, nsegment_table,
     *            block_map_table, nblock_map_table,
     *            address_table, op)
         endif
      endif

      istabvalue = array_table(c_scalar_index,19)
      scalar_table(istabvalue) = stabvalue
C
C Watson Added:
C
      ieom_tol = array_table(c_scalar_index,20)
      scalar_table(ieom_tol) = eom_tol

      ieom_roots = array_table(c_scalar_index,21)
      scalar_table(ieom_roots) = eom_roots

      iexcite = array_table(c_scalar_index,22)
      scalar_table(iexcite) = excite
C
C Watson
C
      next = 0
      next_server = 0
      next_worker = 0
      do i = 23, narray_table   ! 1st 19 arrays are pre-defined.
         type = array_table(c_array_type,i) ! Watson changed from 20 to 22
         if (type .eq. static_array) then

c---------------------------------------------------------------------------
c   Make sure the array is actually being used.
c---------------------------------------------------------------------------

            i_used_it = .false.
            do j = 1, noptable
               if (optable(c_opcode,j) .lt. call_op .or.
     *             optable(c_opcode,j) .gt. sp_ldindex_op) then
                  if (optable(c_op1_array,j) .eq. i) i_used_it = .true.
                  if (optable(c_op2_array,j) .eq. i) i_used_it = .true.
                  if (optable(c_result_array,j) .eq. i) 
     *                      i_used_it = .true.
                  if (i_used_it) go to 1200
               endif
            enddo

 1200       continue
            if (i_used_it) then
               n = 1
               do j = 1, array_table(c_nindex,i)
                  index = array_table(c_index_array1+j-1,i)
                  bseg = index_table(c_bseg, index)
                  eseg = index_table(c_eseg, index)
                  call get_index_segment(index, bseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
                  array_table(c_index_range1+j-1,i) = val1
                  call get_index_segment(index, eseg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
                  array_table(c_index_range2+j-1,i) = val2
                  k = array_table(c_index_range2+j-1,i) -
     *                array_table(c_index_range1+j-1,i) + 1
                  n = n * k
               enddo

               call mem_alloc(xxx, n, size, ixc, .true., ierr)
               if (ierr .ne. 0) then
                  print *,'Mem failure allocating static array #',i
                  call abort_job()
               endif

               address_table(i) = c_loc64(xxx, ixc, size)
               nbytes_static = nbytes_static + n * size

c---------------------------------------------------------------------------
c   Zero out the memory for static arrays.
c---------------------------------------------------------------------------

               do j = 1, n
                  xxx(ixc+j-1) = 0.
               enddo
            endif   ! i_used_it
         endif 

         if (type .eq. local_array) then
c---------------------------------------------------------------------------
c   Save original array indices for use in later block lookups.
c---------------------------------------------------------------------------

            nind = array_table(c_nindex,i)
            do j = 1, nind
               array_table(c_index_original+j-1,i) =
     *               array_table(c_index_array1+j-1,i)
            enddo
         endif

         if (type .eq. distributed_array .or.
     *       type .eq. served_array) then

c-------------------------------------------------------------------------
c   Define the block_map_entry table for the array.
c-------------------------------------------------------------------------

            iproc = 0
            array_table(c_block_map,i) = next+1
            nind = array_table(c_nindex,i)

            nblks = 1
            do k = 1, nind
               ind_arr(k) = array_table(c_index_array1+k-1,i)
               bseg_arr(k) = index_table(c_bseg, ind_arr(k))
               eseg_arr(k) = index_table(c_eseg, ind_arr(k))
               nseg_arr(k) = eseg_arr(k) - bseg_arr(k) + 1 
               nblks = nblks * nseg_arr(k)
            enddo

            do k = nind+1,mx_array_index
               ind_arr(k) = 0
               bseg_arr(k) = 0
               eseg_arr(k) = 0
               nseg_arr(k) = 0
            enddo
 
            array_table(c_numblks,i)     = nblks
            array_table(c_create_flag,i) = 0

c---------------------------------------------------------------------------
c   Save original array indices for use in later block lookups.
c---------------------------------------------------------------------------

            do j = 1, nind
               array_table(c_index_original+j-1,i) = 
     *               array_table(c_index_array1+j-1,i)
            enddo

c----------------------------------------------------------------------------
c   Initialize iseg array.
c---------------------------------------------------------------------------

            do k = 1, nind
               iseg(k) = bseg_arr(k)
            enddo
         
c---------------------------------------------------------------------------
c   Fill in block_map_table info for each blocks of the array.
c---------------------------------------------------------------------------

            do k = 1, nblks
               next = block_map_lookup(iseg, nind, i, array_table(1,i),
     *                    index_table, nindex_table)

               if (type .eq. distributed_array) then
                  block_map_table(c_processor,next) =
     *                               mod(next_worker,ncworkers)
                  next_worker = next_worker + 1
               else
                  next_server = next_server + 1
                  iproc_company_rank = mod(next_server-1,niocompany)
                  block_map_table(c_processor,next) =
     *                           ioproc(iproc_company_rank+1)
               endif

               block_map_table(c_bmap_blkndx,next) = 0
               do j = 1, nind
                  block_map_table(c_block_map_seg+j-1,next) = iseg(j)
               enddo

               do j = nind+1,mx_array_index
                  block_map_table(c_block_map_seg+j-1,next) = 0
               enddo

               j = 1
  100          continue

c-------------------------------------------------------------------------
c   Increment the segment counter of the jth index.
c-------------------------------------------------------------------------

               iseg(j) = iseg(j) + 1
               if (iseg(j) .gt. eseg_arr(j)) then
                  iseg(j) = bseg_arr(j)   ! reset to beginning of index j.

c---------------------------------------------------------------------------
c   Bump the index, repeat if necessary.
c---------------------------------------------------------------------------

                  j = j + 1
                  if (j .le. nind) go to 100
               endif
            enddo   ! k
         endif 

         if (array_table(c_array_type,i) .eq. served_array) then
            if (io_company_id .eq. 0) then
               if (me .eq. 0) 
     *            print *,'Warning: SIAL program uses served arrays,',
     *           ' an IOCOMPANY may be required.'
            endif
         endif 
      enddo

c--------------------------------------------------------------------------
c   Initialize each PARDO instruction with the max. number of batches.
c--------------------------------------------------------------------------

      if (.not. dryrun) then
         call init_stack(noptable)  ! initialize instruction stack. 
         iwhere = 0
         nwhere = 0

         do iop = 1, noptable
            if (optable(c_opcode, iop) .eq. pardo_op) then
               do i = 1, mx_array_index
                  ind_arr(i) = optable(c_ind1+i-1,iop)
               enddo

c---------------------------------------------------------------------------
c   Save where conditions, get any new ones for this pardo.
c---------------------------------------------------------------------------

               ierr = push_do(iwhere, nwhere)
               call get_where_conditionals(optable, iop, noptable)

               max_batch = 0
               call  lb_set_loop_incr_mapping(ind_arr, 
     *             index_table, nindex_table,
     *             segment_table, nsegment_table, 
     *             ind_pardo, bseg_arr, eseg_arr, nind, max_batch)
               optable(c_pardo_max_batch,iop) = max_batch
               
c---------------------------------------------------------------------------
c   Restore previous where conditions.
c---------------------------------------------------------------------------

            else if (optable(c_opcode,iop) .eq. do_op) then
               ierr = push_do(iwhere, nwhere) 
               call get_where_conditionals(optable, iop, noptable) 
            else if (optable(c_opcode,iop) .eq. enddo_op .or.
     *               optable(c_opcode,iop) .eq. endpardo_op) then
               ierr = pop_do(iwhere, nwhere)
            endif
         enddo
      endif   ! .not. dryrun

      if (my_company_rank .eq. 0) then
         print *,'Allocated ',nbytes_static,' bytes for static arrays.'
      endif

      if (my_company_rank .eq. 0) then
c         print *,'INDEX_TABLE:'
c         do i = 1, nindex_table
c            print *,i,': ',(index_table(j,i),j=1,lindex_table_entry)
c         enddo

c         print *,'ARRAY TABLE:'
c         do i = 1, narray_table
c            print *,i,': ',(array_table(j,i),j=1,larray_table_entry)
c         enddo

c         print *,'BLOCK_MAP_TABLE: nblock_map_table ',
c     *           nblock_map_table
c         do i = 1, nblock_map_table
c            print 10101,i,(block_map_table(j,i),j=1,lblock_map_entry)
c         enddo
10101 format(' ',i6,': ',10(i6,1x))
 
c         print *,'OPTABLE:'
c         do i = 1, noptable
c            print *,i,': ',(optable(j,i),j=1,loptable_entry)
c         enddo
 
c         print *,'SEGMENT_TABLE:'
c         do i = 1, nsegment_table
c            print *,i,': ',(segment_table(j,i),j=1,lsegment_table_entry)
c         enddo

c         print *,'ADDRESS TABLE: '
c         do i = 1,narray_table 
c            if (address_table(i) .ne. 0) 
c     *          print *,'Entry ',i,': ',address_table(i)
c         enddo

c         print *,'Initial scalar table:'
c         do i = 1, nscalar_table
c            print *,'Scalar entry ',i,': ',scalar_table(i)
c         enddo
      endif

        call build_proctab(optable, noptable, proctab, nproctab)
	call init_stack(noptable)  ! initialize instruction stack.
      
c-------------------------------------------------------------------------
c   Set up a two timers for each pardo/endpardo:
c   The first is the timer for the pardo loop itself,
c   the second times the amount of block wait time in the pardo loop.
c------------------------------------------------------------------------

      do i = 1, noptable
         optable(c_instr_timer,i) = 0
      enddo

      stat_key = 0
      do i = 1, noptable
         if (optable(c_opcode,i) .eq. pardo_op) then
            timerdesc = ' '
            write (timerdesc,9200) optable(c_lineno,i)
            call register_timer(timerdesc, elapsed_time_timer, key1)

            timerdesc = ' '
            write (timerdesc,9300) optable(c_lineno,i)
            call register_timer(timerdesc, elapsed_time_timer, key2)

            call pack_pardo_timer(optable(c_instr_timer,i),key1,key2) 

c-------------------------------------------------------------------------
c  Set up a server stat_key.
c-------------------------------------------------------------------------

            call translate_key(optable(c_lineno,i), stat_key)
         endif

c-------------------------------------------------------------------------
c   Store the stat_key in the optable.
c   This key is used to sum the server statistics for the pardo.
c-------------------------------------------------------------------------

         optable(c_server_stat_key,i) = stat_key
      enddo

c-------------------------------------------------------------------------
c   Register the timers for each instruction.
c-------------------------------------------------------------------------

      do i = 1, noptable
         opcode = optable(c_opcode,i)

c-------------------------------------------------------------------------
c   Screen out certain opcodes that we don't want to track.  
c   1. assignment (with result type = scalar).
c-------------------------------------------------------------------------

         doit = .true.
         if (opcode .eq. assignment_op) then
            result = optable(c_result_array,i)
            if (array_table(c_array_type,result) .eq. 
     *          scalar_value) doit = .false.
         endif

         if (doit .and. opcode .ne. pardo_op) then
            opdesc = get_operation_desc(opcode)
            if (opdesc(1:1) .ne. ' ') then
               timerdesc = ' '
               write (timerdesc,9000) optable(c_lineno,i),opdesc
               call register_timer(timerdesc, elapsed_time_timer, key)
               optable(c_instr_timer,i) = key
            endif
         endif
      enddo

c-------------------------------------------------------------------------
c   Register a timer for each procedure.  This timer is stored in the
c   "op1" field of each call instruction.
c-------------------------------------------------------------------------

      do i = 1, noptable
         if (optable(c_opcode,i) .eq. call_op) then
            timerdesc = ' '
            write(timerdesc,9100) optable(c_lineno, i)
            call register_timer(timerdesc, elapsed_time_timer, key)
            optable(c_op1_array,i) = key
         endif
      enddo
 
c-------------------------------------------------------------------------
c   Set the "lock index" field for each pardo instruction.  These will be
c   used in case the load-balancing pardo option is used.
c-------------------------------------------------------------------------

      next = 0
      do i = 1, noptable
         if (optable(c_opcode,i) .eq. pardo_op) then
            next = next + 1
            optable(c_pardo_lock_index,i) = next
         endif
      enddo

c      if (my_company_rank .eq. 0) then
c         print *,'Processing Operation table:'
c         do i = 1, noptable
c            print *,'Op #',i,': ',(optable(j,i),j=1,loptable_entry)
c         enddo
c      endif

 9000 format('Line',i6,': ',a12)
 9100 format('Proc call at line ',i6)
 9200 format('Pardo at line ',i6)
 9300 format('Blkwait for pardo ',i6)
      return
      end

