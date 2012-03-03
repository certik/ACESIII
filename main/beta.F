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
      program aces3
c      subroutine beta()
      implicit none

      include 'mpif.h'
      include 'aces3_version.h'
      include 'int_gen_parms.h'
      include 'proto_defines.h'
      include 'proto_events.h'
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'hess.h'
      include 'server_monitor.h'
      include 'timerz.h'
      include 'sial_config_params.h'
      include 'checkpoint_data.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer ierr, mpierr, istat
      integer mstr_is_wrkr
      integer id, platoon_id, nwrkr, mem_wrkr
      integer instance, sial_program_instance
      integer next
      integer i, j, x, y, z, n, dummy
      integer n1, n2
      integer proto_driver
      integer str_trimlen
      integer rank, company, pst_get_company
      integer my_company
      integer pst_get_company_rank
      integer company_comm, pst_get_company_comm
      integer master, pst_get_master
      integer nworkers_needed, nworkers_min
      integer check_malloc
      integer sheap_flag
      integer auto_seg_switch
      integer status(MPI_STATUS_SIZE)
      integer company_id
      integer saveocc, savevirt, saveao

      character*20 logfile
      character*80 sial_program
      character*80 company_string
      character*80 package
      character*80 tracetype
      character*80 erd_fast_memcalc
      character*80 debug_print
      character*80 restart_program
      character*80 auto_seg_algorithm
      character*80 timer_param
      character*80 ctrip_str
      character*80 ihess_str
      character*80 jhess_str
      character*80 sub_str

      logical company_defined
      logical is_file_open
      logical*8 l8false
      integer jodadone
      logical first
      logical instab
      logical really_ignore_dropmo
      logical scfskip
      logical dryrun_abort
      logical do_restart, file_exist
      logical flag

      double precision t1, t2
      double precision vnn, diff, avg
      double precision dbuf(1)

      integer i1, i2
      integer*8 mpierr8, nprocs8

#ifdef ALTIX
      pointer (dptr, dbuf)
#else
      common dbuf
#endif

      call init_machine_types()
      call aces_init_rte()
      call mpi_init(ierr)

#ifdef ALTIX
      sheap_flag = 1
#else
      sheap_flag = 0
#endif

      call mpi_comm_size(mpi_comm_world, nprocs, mpierr)
      call mpi_comm_rank(mpi_comm_world, me, mpierr)

c---------------------------------------------------------------------------
c   Scan the ZMAT file for the MAXMEM parameter.  This is needed to set
c   up the heap memory that will be used.
c---------------------------------------------------------------------------

      if (me .eq. 0) then
         call init_params('ZMAT')
         maxmem = 0
         call igetparam('MAXMEM', 1, maxmem)
      endif

c---------------------------------------------------------------------------
c   Initial call to set up memory.
c---------------------------------------------------------------------------

      call mpi_bcast(maxmem, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      if (maxmem .eq. 0) maxmem = check_malloc() * 0.8
      call mem_alloc_init(maxmem, sheap_flag, ierr)
      if (ierr .lt. 0) then
         print *,'MEMORY ALLOCATION ERROR: maxmem = ',maxmem
         call abort_job()
      endif
#ifdef ALTIX
      dptr = dshptr    ! set dbuf pointer to shared mem.
#endif

      dryrun_abort = .false.
      if (me .eq. 0) then
         print *,'ACESIII version, ',aces3_version,
     *            ' is run on ',nprocs,' processors.'
c         call prt_time('Program start')

c-------------------------------------------------------------------------
c   Do we have an existing JOBARC file?  If not, we are doing an SCF
c   from scratch.
c-------------------------------------------------------------------------

         call scf_init()

c-------------------------------------------------------------------------
c   Verify that we have valid combinations of parameters from the *ACES2
c   section of the ZMAT file.
c-------------------------------------------------------------------------

        call verify_aces2_params()

c-------------------------------------------------------------------------
c   Create a ZMAT.AUTO file from the user's ZMAT, adding the defaulted
c   parameters in the *SIP section.
c-------------------------------------------------------------------------

         call write_zmat_auto()
         call aces_ja_fin()   ! close JOBARC subsystem

c--------------------------------------------------------------------------
c   Initialize the parameter decoding routines.
c--------------------------------------------------------------------------

         call init_params('ZMAT.AUTO')

c-------------------------------------------------------------------------
c   Scan all sio files to determine whether second derivative integrals
c   are needed.
c-------------------------------------------------------------------------
         stabvalue = 0.0d0
         instab = .true.
         scfskip = .false.
         really_ignore_dropmo = .false.
         sial_program_instance = 1
         first = .true.
   10    continue
         sial_program = ' '
         call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                  sial_program)
         if (sial_program .ne. ' ') then
            call get_sial_config_params(sial_program)
            if (.not. use_2der_integrals) then
               sial_program_instance = sial_program_instance + 1
	       go to 10
            endif
         endif
         print *,'USE_2DER_INTEGRALS = ',use_2der_integrals
      endif

      call mpi_barrier(mpi_comm_world, mpierr)

      calc_2der_integrals = use_2der_integrals

      sial_program_instance = 0
      jodadone = 0
      dryrun   = .true.   ! 1st pass through SIAL code is a dryrun
      simulator = .false.
      nworkers_min = 0
      server_monitor_on = .false.

c-----------------------------------------------------------------------
c   Check for existence of the restart file.
c-----------------------------------------------------------------------

      restart_file        = 94
      inquire (file='aces3.restart', exist = do_restart)

c------------------------------------------------------------------------
c   Initialize checkpoint_data common block.
c------------------------------------------------------------------------

      master_ckpt_unit    = 93
      nckpt_arrays        = 0
      nactive_allocate_table = 0
      nactive_create_table   = 0
      nfree_space            = 0
      nfree_space_candidates = 0
      restart_job = do_restart

      call mpi_barrier(mpi_comm_world, mpierr)

   50 continue
      call reset_internal_system_data()   ! restore system flags
      write_hess = .false.
      t1 =  mpi_wtime()

      if (me .eq. 0) then

c---------------------------------------------------------------------
c   Set parameter defaults.
c---------------------------------------------------------------------

         mstr_is_wrkr   = 1
         package        = 'ERD'
         molfile     = 'MOL'
         local_path     = '.'
         tracetype      = ' '
         sial_program   = ' '
         erd_fast_memcalc    = 'NO'
         debug_print         = 'NO'
         restart_program = ' '
         auto_seg_algorithm = 'SEGMENT_OPTIMIZED'
         timer_param        = 'NO'
         aat_database = ' '
         ctrip_str    = ' '
         ihess_str    = ' '
         jhess_str    = ' '
         sub_str      = ' '

         nworkthread = 20
         sip_mx_ao_segsize   = 0
         sip_mx_occ_segsize  = 0
         sip_mx_virt_segsize = 0
         sip_sub_segsize     = 5
         sip_sub_occ_segsize = 5
         sip_sub_virt_segsize = 5
         sip_sub_ao_segsize   = 5
         reg = 0.0d0

c---------------------------------------------------------------------
c   Read parameters.
c---------------------------------------------------------------------

         call igetparam('MASTER_IS_WORKER', 1, mstr_is_wrkr)
         call igetparam('SIP_MX_SEGSIZE', 1, sip_mx_ao_segsize)
         call igetparam('SIP_MX_OCC_SEGSIZE', 1, sip_mx_occ_segsize)
         call igetparam('SIP_MX_VIRT_SEGSIZE', 1, sip_mx_virt_segsize)
         call cgetparam('INTEGRAL_PACKAGE', 1, package)
         call cgetparam('MOLFILE', 1, molfile)
         call cgetparam('LOCAL_PATH', 1, local_path)
         call cgetparam('TRACE', 1, tracetype)
         call cgetparam('FAST_ERD_MEMCALC', 1, erd_fast_memcalc)
         call cgetparam('DEBUG_PRINT', 1, debug_print)
         call igetparam('NWORKTHREAD', 1, nworkthread)
         call cgetparam('AUTO_SEG_ALGORITHM', 1, auto_seg_algorithm)
         call cgetparam('AAT_DATABASE', 1, aat_database)
         call cgetparam('ITRIP', 1, ctrip_str)
         call cgetparam('IHESS', 1, ihess_str)
         call cgetparam('JHESS', 1, jhess_str)
         call cgetparam('SUB', 1, sub_str)
         call igetparam('SIP_SUB_SEGSIZE', 1, sip_sub_segsize)
         call igetparam('SIP_SUB_OCC_SEGSIZE', 1, sip_sub_occ_segsize)
         call igetparam('SIP_SUB_VIRT_SEGSIZE', 1, sip_sub_virt_segsize)
         call igetparam('SIP_SUB_AO_SEGSIZE', 1, sip_sub_ao_segsize)
         call rgetparam('REGULARIZER',1,reg)


         call decode_fmo_params()

c-------------------------------------------------------------------------
c   Strip any embedded nulls out of molfile name.
c-------------------------------------------------------------------------

         do i = 1, len(molfile)
            if (molfile(i:i) .eq. char(0)) molfile(i:i) = ' '
         enddo

         ierr = 0
         if (debug_print(1:1) .eq. 'Y') then
            dbg = .true.
         else if (debug_print(1:1) .eq. 'N') then
           dbg = .false.
         else
           print *,'====>> Invalid value for DEBUG_PRINT'
           ierr = ierr + 1
         endif

c-----------------------------------------------------------------------
c   Open the restart file.
c-----------------------------------------------------------------------

         open (unit=restart_file, file = 'aces3.restart')

         if (first) then
            maxmem = 0
            call igetparam('MAXMEM', 1, maxmem)
            if (do_restart) then
               read (restart_file, *) restart_program

               if (restart_program(1:1) .ne. ' ') then
                  n1 = str_trimlen(restart_program)
                  print *,'***** Restarting job on ',
     *                 restart_program(1:n1),'*****'
                  dryrun     = .false.    ! skip dryrun pass on restart

c---------------------------------------------------------------------------
c   Find the restart_program in the list of SIAL_PROGRAMs.
c---------------------------------------------------------------------------

                  do i = 1, 100
                     sial_program = ' '
                     call cgetparam('SIAL_PROGRAM', i,
     *                     sial_program)
                     if (sial_program(1:1) .eq. ' ') then
                        print *,'Error: RESTART program not found in ',
     *                    'list of SIAL_PROGRAM parameters.'
                        call abort_job()
                     else
                        n2 = str_trimlen(sial_program)
                        if (sial_program(1:n2) .eq.
     *                      restart_program(1:n1)) then
                           sial_program_instance = i - 1 ! step back one program
                           go to 60
                        endif
                     endif
                  enddo
  60              continue
               endif    ! restart_program
             endif      ! do_restart
         endif          ! first

c---------------------------------------------------------------------------
c   Get the next instance of the SIAL_PROGRAM parameter.
c---------------------------------------------------------------------------
         if (instab .and. (.not. dryrun) ) then
c Check to see if there is an instability to follow
C Check to see if we just ran the instability sial program
           sial_program = ' '
           call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                    sial_program)
           if ((sial_program(1:6) .eq. 'instab') .or.
     *         (sial_program(1:6) .eq. 'INSTAB')) then
             sial_program_instance = 1
C If we just ran instab and there was no negative eigenvalue,
C then stop checking
             if (stabvalue .ge. 0.0d0) instab = .false.
             scfskip = .true.
           endif
C turn off the really_ignore_dropmo flag before we check to see if we need it
           really_ignore_dropmo = .false.
         endif
         sial_program_instance = sial_program_instance + 1
         sial_program = ' '
         call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                  sial_program)
         if (instab) then
           if ((sial_program(1:4) .eq. 'tran') .or.
     *         (sial_program(1:4) .eq. 'TRAN')) then
C If the next program is an integral transform and then AFTER that an instability
C check, then make sure NOT to drop mos (turn really_ignore_dropmo = .true.)
               sial_program = ' '
               call cgetparam('SIAL_PROGRAM', sial_program_instance+1,
     *                        sial_program)
             if ((sial_program(1:6) .eq. 'instab') .or.
     *           (sial_program(1:6) .eq. 'INSTAB'))
     *           really_ignore_dropmo = .true.
             sial_program = ' '
             call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                      sial_program)
           endif
         endif
         if ((.not. instab).and. (.not. dryrun) ) then
C If we have stopped checking for instabilities and the instab program is
C called, ignore it
           if ((sial_program(1:6) .eq. 'instab') .or.
     *         (sial_program(1:6) .eq. 'INSTAB')) then
             sial_program_instance = sial_program_instance + 1
             sial_program = ' '
             call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                      sial_program)
           endif
         endif
         if (sial_program(1:1) .eq. ' ') then
            if (dryrun) then
               if (dryrun_abort) then
                  print *,'Job requires at least ',nworkers_min,
     *              ' processors in COMPANY parameter to run.'
                  print *,'----- Aborting due to dryrun errors. ----'
                  jodadone = 1
               endif

               dryrun = .false.
               scfskip = .false.
            else

c---------------------------------------------------------------------------
c   We have looped over the number of SIAL_PROGRAM instances.  Now we
c   must let xjoda do its processing, checking the JODADONE flag to
c   determine whether the job is complete.
c---------------------------------------------------------------------------

               if (vib_freq_calc .or. geom_opt) then
                  call mem_alloc_reset()            ! geopt may use some heap
                  call check_jodadone(jodadone)
                  if (sial_program_instance .eq. 1) then
                     print *,'Error: No SIAL_PROGRAM parameter'
                     call abort_job()
                  endif
                  if (jodadone .eq. 0) instab = .true.
               else
                  jodadone = 1
               endif   ! vib_freq

            endif   ! dryrun
         endif   ! sial_program
      endif    ! me .eq. 0

c----------------------------------------------------------------------------
c   Send jodadone flag to all processors, check for end of job.
c---------------------------------------------------------------------------

      call c_flush_stdout()  ! flush joda print to stdout
      call mpi_bcast(jodadone, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      call mpi_barrier(mpi_comm_world,mpierr)
      if (jodadone .ne. 0) then
         if (me .eq. 0) then
            call prt_time('Program stop')

c----------------------------------------------------------------------------
c   Remove the restart file.
c----------------------------------------------------------------------------

            close (restart_file)
            call f_deletefile('aces3.restart' // char(0), -1)
         endif

         call mpi_finalize(ierr)
         stop
c         return
      endif

c--------------------------------------------------------------------------
c   Set up memory allocation subsystem.
c--------------------------------------------------------------------------

      call mpi_bcast(maxmem, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      if (maxmem .eq. 0) maxmem = check_malloc() * 0.8
c      call mem_alloc_init(maxmem, sheap_flag, ierr)
      if (ierr .lt. 0) then
         print *,'MEMORY ALLOCATION ERROR: maxmem = ',maxmem
         call abort_job()
      endif
#ifdef ALTIX
      dptr = dshptr    ! set dbuf pointer to shared mem.
#endif

c--------------------------------------------------------------------------
c   Resume parameter checking on process 0.
c--------------------------------------------------------------------------

      if (me .eq. 0) then
         if (sial_program(1:1) .eq. ' ') then
            sial_program_instance = 1
            call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                  sial_program)
            if ((sial_program(1:3) .eq. 'scf') .or.
     *          (sial_program(1:3) .eq. 'SCF')) then
                if (scfskip) then
                  sial_program_instance = sial_program_instance + 1
                  call cgetparam('SIAL_PROGRAM', sial_program_instance,
     *                           sial_program)
                endif
            endif
         endif

         if (sial_program(1:1) .eq. ' ') then
            print *,'====>> SIAL_PROGRAM parameter is required.'
            ierr = ierr + 1
         else
            print *,'Next SIAL program: ',sial_program
         endif

c--------------------------------------------------------------------------
c   Read program-specific parameters from the config file.
c--------------------------------------------------------------------------

         call get_sial_config_params(sial_program)
         if (really_ignore_dropmo) ignore_dropmo = .true.

c------------------------------------------------------------------------
c   Write current sial_program into the restart file.
c------------------------------------------------------------------------

         rewind (restart_file)
         write (restart_file,*) sial_program
         close (restart_file)

         if (mstr_is_wrkr .eq. 0) then
            master_is_worker = .false.
         else
            master_is_worker = .true.
         endif

         managers_are_workers = .true.

c----------------------------------------------------------------------
c   Set up companies.
c----------------------------------------------------------------------

         call init_company
         company_defined = .false.
         instance = 1

  100   continue
        company_string = ' '
        call cgetparam('COMPANY', instance, company_string)

        if (company_string(1:1) .ne. ' ') then
           read (company_string, *) id, platoon_id, nwrkr,
     *                               mem_wrkr
           if (mem_wrkr .eq. 0) mem_wrkr = maxmem

           if (no_servers) then
              print *,'Overriding COMPANY parameter, using ',
     *            nprocs,' processes as workers.'
              nwrkr = nprocs
           endif

           call build_company_platoon(id, platoon_id,
     *                      nwrkr, mem_wrkr, .false., sial_program)
           if (dbg) print *,'Built company ',id,' platoon ',platoon_id,
     *                       '...'
           instance = instance + 1
           company_defined = .true.
           go to 100
        endif

        if (.not. company_defined) then
           print *,'====>> COMPANY parameter was not specified.'
           ierr = ierr + 1
        endif

        company_string = ' '
        io_company_id = 0

        if (.not. no_servers) then
           call cgetparam('IOCOMPANY', 1, company_string)
           n = str_trimlen(company_string)
           if (n .gt. 0) then
              read (company_string, *) id, platoon_id, nwrkr,
     *                               mem_wrkr
              if (mem_wrkr .eq. 0) mem_wrkr = maxmem
              call build_company_platoon(id, platoon_id,
     *                      nwrkr, mem_wrkr, .true., ' ')
              io_company_id = id
           endif
        endif

        if (package(1:3) .eq. 'ERD') then
           intpkg = flocke_package
        else
           print *,'====>> Invalid integral package: ',package
           ierr = ierr + 1
        endif

        if (erd_fast_memcalc(1:1) .eq. 'Y') then
           fast_erd_memcalc = .true.
        else if (erd_fast_memcalc(1:1) .eq. 'N') then
           fast_erd_memcalc = .false.
        else
           print *,'====>> Invalid value for FAST_ERD_MEMCALC'
           ierr = ierr + 1
        endif

        if (auto_seg_algorithm(1:16) .eq. 'MEMORY_OPTIMIZED') then
           auto_seg_switch = 1
        else if (auto_seg_algorithm(1:17) .eq.
     *                                     'SEGMENT_OPTIMIZED') then
           auto_seg_switch = 2
        else
           print *,'====>> Invalid value for AUTO_SEG_ALGORITHM'
           ierr = ierr + 1
        endif

c----------------------------------------------------------------------------
c   Check timer parameter.  The values can be either 'YES', 'NO' or a
c   sial program.
c----------------------------------------------------------------------------

        next = 0
        n = str_trimlen(sial_program)
  200   continue
        timer_param = ' '
        next = next + 1
        call cgetparam('TIMERS', next, timer_param)
        n2 = str_trimlen(timer_param)

        do_timer = .false.
        if (timer_param(1:1) .eq. 'Y') then
           do_timer = .true.
        else if (timer_param(1:1) .ne. ' ') then

c---------------------------------------------------------------------------
c   Does the timer_param match the sial_program?
c---------------------------------------------------------------------------

           if (timer_param(1:n2) .eq. sial_program(1:n)) then
              do_timer = .true.
           else
              go to 200   ! check next possible value.
           endif
        endif

        tracelevel = 0
        if (tracetype(1:12) .eq. 'CONTRACTIONS') then
           tracelevel = contraction_trace
        else if (tracetype(1:12) .eq. 'INSTRUCTIONS') then
           tracelevel = instruction_trace
        else if (tracetype(1:5) .eq. 'PROCS') then
           tracelevel = proc_trace
        endif

        iworker_type = direct_integral_worker
        cworker_type = sip_contraction_worker
        algorithm_flag = 3

c-------------------------------------------------------------------
c   Read basis data from the molfile.
c-------------------------------------------------------------------

         if (sial_program(1:3) .eq. 'scf' .or.
     *       sial_program(1:3) .eq. 'SCF') then
            if_scf = .true.
#ifdef AIX
            stack_algorithm_type = 1   ! avoids comm subsystem mem exhaustion
#else
            stack_algorithm_type = 2
#endif
         else
            if_scf = .false.
            stack_algorithm_type = 1
         endif

         if (sial_program(1:5) .eq. 'guess' .or.
     &        sial_program(1:5) .eq. 'GUESS') then
           init_scf = .true.
#ifdef AIX
            stack_algorithm_type = 1   ! avoids comm subsystem mem exhaustion
#else
            stack_algorithm_type = 2
#endif
         else
           init_scf = .false.
           stack_algorithm_type = 1
         endif

         first = .false.
         call init_data()
         call joda_prep()

         if (sip_mx_occ_segsize.lt.1.or.
     &       sip_mx_virt_segsize.lt.1.or.
     &       sip_mx_ao_segsize.lt.1) then
c         o This will crash if it cannot find segment sizes.
           saveocc  = sip_mx_occ_segsize
           savevirt = sip_mx_virt_segsize
           saveao   = sip_mx_ao_segsize
           call det_segsizes(auto_seg_switch, nprocs)
           if (saveocc .gt. 0) sip_mx_occ_segsize = saveocc
           if (savevirt .gt. 0) sip_mx_virt_segsize = savevirt
           if (saveao .gt. 0) sip_mx_ao_segsize = saveao
         end if

c----------------------------------------------------------------------
c   Define some SIAL-accessible parameter values.
c----------------------------------------------------------------------

         if (ctrip_str .eq. ' ') then
            itrips = 1
            itripe = max(nalpha_occupied, nbeta_occupied)
         else
            read (ctrip_str, *) i1, i2
            itrips = i1
            itripe = i2
         endif

         if (ihess_str .eq. ' ') then
            ihess1 = 1
            ihess2 = ncenters
         else
            read (ihess_str, *) i1, i2
            ihess1 = i1
            ihess2 = i2
         endif

         if (jhess_str .eq. ' ') then
            jhess1 = 1
            jhess2 = ncenters
         else
            read (jhess_str, *) i1, i2
            jhess1 = i1
            jhess2 = i2
         endif

         if (sub_str .eq. ' ') then
            subb = -1
            sube = -1
         else
            read (sub_str, *) i1, i2
            subb = i1
            sube = i2
         endif

c----------------------------------------------------------------------
c   Additional error checking on parameter values.
c----------------------------------------------------------------------

        if (ncenters .gt. max_centers) then
           print *,'Number of centers is ',ncenters,', greater than ',
     *             ' the limit of ',max_centers
           ierr = ierr + 1
        endif

        if (nshells .gt. max_shells) then
           print *,'Number of shells is ',nshells,', greater than ',
     *             ' the limit of ',max_shells
           ierr = ierr + 1
        endif

        if (ierr .gt. 0) then
           print *,'Task ',me,
     *           ' ====>> Aborting due to parameter errors...'
           call c_flush_stdout()
           call mpi_abort(mpi_comm_world, 1, mpierr)
        else

           call int_gen_init()

           print *,'----------- Parameter Summary ----------------'
           print *,'DRYRUN...........................',dryrun
           print *,'INTEGRAL_PACKAGE.................',package
           print *,'MOLFILE..........................',molfile
           print *,'MAXMEM...........................',maxmem
           print *,'OCCUPIED.........................',noccupied
           print *,'ALPHA_OCCUPIED...................',nalpha_occupied
           print *,'BETA_OCCUPIED....................',nbeta_occupied
           print *,'AO SEGMENT SIZE..................',
     *                     sip_mx_ao_segsize
           print *,'OCCUPIED SEGMENT SIZE............',
     *                     sip_mx_occ_segsize
           print *,'VIRTUAL SEGMENT SIZE.............',
     *                     sip_mx_virt_segsize
           print *,'Local scratch directory..........',local_path
           print *,'Number of centers................',ncenters
           print *,'Number of shells.................',nshells
           print *,'Number of basis functions........',
     *                 end_nfps(nshells)
           print *,'NO_SERVERS.......................',no_servers
           print *,'IGNORE DROPMO PARAMETERS.........',
     *                                        ignore_dropmo
           print *,'ITRIP............................',itrips,itripe
           print *,'IHESS............................',ihess1,ihess2
           print *,'JHESS............................',jhess1,jhess2
           print *,'SIP_SUB_SEGSIZE..................',sip_sub_segsize
           print *,'SIP_SUB_OCC_SEGSIZE..............',
     *                                            sip_sub_occ_segsize
           print *,'SIP_SUB_VIRT_SEGSIZE.............',
     *                                            sip_sub_virt_segsize
           print *,'SIP_SUB_AO_SEGSIZE...............',
     *                                            sip_sub_ao_segsize
           print *,'REGULARIZER......................',reg
           if (subb .ne. -1 .or. sube .ne. -1)
     *        print *,'SUB..............................',subb,sube
           if (do_timer)
     *        print *,'An instruction timer report will be printed'
           print *,' '
           if (dbg) then
              print *,'end_nfps: ',(end_nfps(i),i=1,nshells)
           endif

        endif
      else
         call aces_init_rte()     ! initialize some aces common blocks
      endif

c----------------------------------------------------------------------
c   Wait for rank 0 before proceeding.
c----------------------------------------------------------------------

      call mpi_barrier(mpi_comm_world, mpierr)

c---------------------------------------------------------------------
c   Distribute the common block to all processors.
c---------------------------------------------------------------------

      call int_send_common()

c---------------------------------------------------------------------
c   Make sure the tracelevel is set on all procs.
c---------------------------------------------------------------------

      call mpi_bcast(tracelevel, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      call mpi_bcast(dryrun, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      call mpi_bcast(dbg, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      call mpi_bcast(do_timer, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)
      call mpi_bcast(vvvi_stack, 1, mpi_integer, 0,
     *               mpi_comm_world, mpierr)

c--------------------------------------------------------------------
c   Do the event processing.
c--------------------------------------------------------------------

      call build_pst(nprocs, master_is_worker, managers_are_workers,
     *               maxmem, io_company_id)
      if (me .eq. 0) then

c--------------------------------------------------------------------------
c   Send the basis function data to each worker.
c   Send the transformation coeffs and eigenvalues to each non-I/O company
c   worker.
c--------------------------------------------------------------------------

         my_company = pst_get_company(me)
         company_comm = pst_get_company_comm(me)
         my_company_rank = pst_get_company_rank(me)

         do i = 1, nprocs
            company      = pst_get_company(i-1)
            rank = pst_get_company_rank(i-1)
            if (company .eq. my_company .and.
     *          rank .ne. my_company_rank) then
               call mpi_isend(dbuf(master_icoord), nalloc_pkg,
     *                  mpi_double_precision, rank,
     *                  basis_function_request_event,
     *                  company_comm, dummy, mpierr)
               call enter_request(dummy)
            endif
         enddo
      else

c------------------------------------------------------------------------
c   Zero the gradient and hessian on all except the master.
c------------------------------------------------------------------------

         do i = 1, 3*ncenters
            gradient_data(i) = 0.
         enddo


         call mem_alloc(dbuf, 9*ncenters*ncenters,
     *               bytes_per_double, ix_hess,
     *               .true., ierr)
         if (ierr .ne. 0) then
            print *,'Error attempting to malloc hessian data.'
            call mpi_abort(mpi_comm_world, 1, ierr)
         endif

         do i = 1, 9*ncenters*ncenters
            dbuf(ix_hess+i-1)   = 0.d0
         enddo
      endif

c--------------------------------------------------------------------------
c   Compile the SIAL program(s) on the master processor, distribute them
c   to the appropriate companies.
c--------------------------------------------------------------------------

#ifdef HP
      call sial_compile_and_send_tables()
#else
      call sial_compile_and_send_tables1()
#endif

      call mpi_barrier(mpi_comm_world, mpierr)
      call worker_init(istat, nworkers_needed)
      if (istat .ne. 0)
     *   nworkers_min = max(nworkers_min, nworkers_needed)
      call prt_time('Returned from worker_init')
      call worker_work()
      call prt_time('Returned from worker_work')
      call worker_termination()
      call prt_time('Returned from worker_termination')

c--------------------------------------------------------------------------
c   Print out the Nuclear-Nuclear gradient and determine the symmetry
c   equivalence of atoms.
c--------------------------------------------------------------------------

      if (me .eq. 0)  then
         if (dryrun) then
            if (istat .eq. 0) then
               print *,'*** COMPLETED DRYRUN PASS FOR ',
     *           sial_program
            else
               print *,'*** FAILED IN DRYRUN FOR STEP ',
     *           sial_program
               dryrun_abort = .true.
            endif
            call aces_ja_fin()
            go to 1000  ! skip this if dryrun
         endif

        call handle_gradient()
        if (write_hess)
     *   call handle_hessian(dbuf(ix_hess), dbuf(ix_nnhess),
     *                       ncenters, write_hess)

c---------------------------------------------------------------------------
c   Write SCFENEG, TOTENERG, GRADIENT, and HESSIAN data to JOBARC.
c   Clear the dirty flag on JOBARC, and close the JOBARC files.
c---------------------------------------------------------------------------

            call dputrec(1,'JOBARC', 'SCFENEG', 1, scf_energy)
            call dputrec(1,'JOBARC', 'TOTENERG', 1, totenerg)
            call dputrec(1, 'JOBARC', 'GRADIENT', 3*ncenters,
     *                   gradient_data)
            call iputrec(1, 'JOBARC', 'NAOBASFN', 1, naobasis)
            call iputrec(1, 'JOBARC', 'NBASTOT', 1, nbasis)
            call iputrec(1, 'JOBARC', 'AIIIJARC', 1, 1)
            if (write_hess) then
               call dputrec(1, 'JOBARC', 'HESSIANM',
     *            (3*ncenters)*(3*ncenters),
     *            dbuf(ix_hess))
               print *,'WRITE HESSIANM TO JOBARC: LENGTH = ',
     *          (3*ncenters)*(3*ncenters)
            endif
            call iputrec(1,'JOBARC','DIRTYFLG',1,0)

            call aces_ja_fin()
            if (dbg) print *, 'CALLED ACES_JA_FIN'

            if (dbg) print *,
     *         '----------- Normal termination of beta -----------'
            t2 = mpi_wtime()
            print *,'Elapsed time = ',t2-t1,' seconds.'

c---------------------------------------------------------------------
c   Receive timer data, print reports.
c---------------------------------------------------------------------

            if (pst_get_company(me) .ne. io_company_id)
     *             call free_all_blocks()
            call mem_alloc_reset()
            call tmanal(sial_program, io_company_id)
      endif

1000  continue
      call mpi_barrier(mpi_comm_world, mpierr)

c--------------------------------------------------------------------------
c   Clear any leftover MPI requests that were made earlier.
c--------------------------------------------------------------------------

c      if (me .eq. 0 .and. .not. dryrun)
c     *      call clear_static_data_requests()

      if (pst_get_company(me) .ne. io_company_id)
     *           call free_all_blocks()
      call mem_alloc_reset()
      call mpi_barrier(mpi_comm_world, mpierr)
      call free_mpi_data()
      call c_flush_stdout()

c----------------------------------------------------------------------------
c   Remove checkpoint files if they exist.
c----------------------------------------------------------------------------

      if (me .eq. 0 .and. .not. dryrun) then
         inquire (file='master_ckpt_file', exist = file_exist)
         if (file_exist) 
     *      call f_deletefile('master_ckpt_file'//char(0),-1)
         inquire (file='CKPT_DIST',exist = file_exist)
         if (file_exist)
     *      call f_deletefile('CKPT_DIST'//char(0),-1)
      endif
      call mpi_barrier(mpi_comm_world, mpierr)
      go to 50

      end

