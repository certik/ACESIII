c--------------------------------------------------------------------------
c   Common block saved_data contains data that must be reinitialized
c   for the execution of a new SIAL program.  
c
c   This reinitialization will be done in subroutine 
c   reset_internal_system_data
c----------------------------------------------------------------------------
      common /saved_data/first_master_send_work,
     *           first_sip_get_next_work,
     *           terminated_worker_termination,
     *           printwarn_worker_work, flag_internal_event,
     *           master_error_processing, master_called,
     *           managers_notified,
     *           first_timer, first_aces_to_erd,
     *           first_contract_blocks, first_dcont2,
     *           first_der_int_setup, first_operation_desc,
     *           first_blk_event_log, cwork_debug,
     *           first_diis_setup,
     *           first_blocks_to_list, first_list_to_blocks,
     *           n_saved_requests, cwork_iprint,
     *           i_been_called_list_to_blocks,
     *           i_been_called_shmem_list

      logical first_master_send_work           ! set to .true.
      logical first_sip_get_next_work          ! set to .true.
      logical terminated_worker_termination    ! set to .false.
      logical printwarn_worker_work            ! set to .true.
      logical flag_internal_event              ! set to .false.
      logical master_error_processing          ! set to .false.
      logical master_called                    ! set to .false.
      logical managers_notified                ! set to .false.
      logical first_timer                      ! set to .true.
      logical first_aces_to_erd                ! set to .true.
      logical first_contract_blocks            ! set to .true.
      logical first_dcont2                     ! set to .true.
      logical first_der_int_setup              ! set to .true.
      logical first_operation_desc             ! set to .true.
      logical first_blk_event_log              ! set to .true.
      logical cwork_debug                      ! set to .false.
      logical first_diis_setup                 ! set to .true.
      logical first_blocks_to_list             ! set to .true.
      logical first_list_to_blocks             ! set to .true.
      integer n_saved_requests                 ! 0
      integer cwork_iprint                     ! 0
      integer i_been_called_list_to_blocks     ! 0
      integer i_been_called_shmem_list         ! 0
