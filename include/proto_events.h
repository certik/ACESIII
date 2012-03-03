
c------------------------------------------------------------------------
c   Definition of prototype events
c------------------------------------------------------------------------

       integer null_event
       integer master_event_tag
       integer worker_event_tag
       integer manager_event_tag
       integer status_request_tag
       integer status_response_tag
       integer status_set_request_tag
       integer timer_data_request_tag
       integer timer_desc_request_tag
       integer integral_worker_identity_tag


       parameter (null_event      = -1)
       parameter (master_event_tag       = 9999)
       parameter (status_request_tag     = 9998)
       parameter (status_response_tag    = 9997)
       parameter (status_set_request_tag = 9996)
       parameter (timer_data_request_tag = 9995)
       parameter (integral_worker_identity_tag = 9994)
       parameter (worker_event_tag       = 9992)
       parameter (manager_event_tag      = 9991)
       parameter (timer_desc_request_tag = 9990)

c--------------------------------------------------------------------------
c   Worker-targeted events.
c--------------------------------------------------------------------------

       integer wrkr_abort_event
       integer wrkr_terminate_event 
       integer wrkr_initialization_event

       parameter (wrkr_abort_event          = 1111)
       parameter (wrkr_terminate_event      = 1112)
       parameter (wrkr_initialization_event = 1113)

c---------------------------------------------------------------------------
c   Master-targeted events.
c---------------------------------------------------------------------------

       integer mstr_initialization_event
       integer mstr_terminate_wrkr_event
       integer mstr_abort_event
       integer mstr_terminate_mgr_event
       integer mstr_init_mgr_response_event
       integer mstr_init_wrkr_response_event
       integer mstr_work_request_event
       integer mstr_error_event

       parameter (mstr_initialization_event = 2222)
       parameter (mstr_terminate_wrkr_event = 2223)
       parameter (mstr_abort_event          = 2224)
       parameter (mstr_terminate_mgr_event = 2225)
       parameter (mstr_init_mgr_response_event  = 2226)
       parameter (mstr_init_wrkr_response_event  = 2227)
       parameter (mstr_work_request_event        = 2228)
       parameter (mstr_error_event               = 2229)

c---------------------------------------------------------------------------
c   Master-targeted events sent by the contraction worker mission.
c---------------------------------------------------------------------------

       integer mstr_integral_request_event
       parameter (mstr_integral_request_event       = 2230)
       integer basis_function_request_event
       parameter (basis_function_request_event      = 2231)
       integer scf_coeff_request_event
       parameter (scf_coeff_request_event           = 2232)
       integer scfa_coeff_request_event
       parameter (scfa_coeff_request_event          = 2233)
       integer scfb_coeff_request_event
       parameter (scfb_coeff_request_event          = 2234)
       integer epsa_coeff_request_event
       parameter (epsa_coeff_request_event          = 2235)
       integer epsb_coeff_request_event
       parameter (epsb_coeff_request_event          = 2236)
       integer predef_constant_tag
       parameter (predef_constant_tag               = 2237)
       integer focka_coeff_request_event
       parameter (focka_coeff_request_event         = 2238)
       integer fockb_coeff_request_event
       parameter (fockb_coeff_request_event         = 2239)

c---------------------------------------------------------------------------
c   Master-targeted events used in communication with the SIP server.
c---------------------------------------------------------------------------

       integer mstr_sip_prepare_event
       parameter (mstr_sip_prepare_event   = 2333)
       integer mstr_sip_request_event
       parameter (mstr_sip_request_event   = 2334)

c---------------------------------------------------------------------------
c   Manager-targeted events.
c---------------------------------------------------------------------------

       integer mgr_abort_event
       integer mgr_terminate_event
       integer mgr_initialization_event

       parameter (mgr_abort_event          = 3553)
       parameter (mgr_terminate_event      = 3554)
       parameter (mgr_initialization_event = 3555)

c-------------------------------------------------------------------------
c   Manager events defined for the "integral mission".
c-------------------------------------------------------------------------

       integer mgr_local_memory_store_event
       integer mgr_local_disk_store_event
       parameter (mgr_local_memory_store_event = 3556)
       parameter (mgr_local_disk_store_event   = 3557)
       
c-------------------------------------------------------------------------
c   Worker events defined for the "integral mission".
c-------------------------------------------------------------------------

      integer wrkr_integral_request_event
      integer wrkr_integral_response_event
      integer wrkr_delay_event
      integer ram_based_work_event
      parameter (wrkr_integral_request_event = 4444)
      parameter (wrkr_integral_response_event = 4445)
      parameter (wrkr_delay_event = 4446)
      parameter (ram_based_work_event = 4443) 

c------------------------------------------------------------------------
c   Worker events defined for the contraction worker mission.
c------------------------------------------------------------------------

      integer wrkr_contraction_request_event
      parameter (wrkr_contraction_request_event  = 4447)
      integer wrkr_contraction_rotation_event 
      parameter (wrkr_contraction_rotation_event  = 4448)

c------------------------------------------------------------------------
c   Worker events defined for the "run_aces" mission.
c------------------------------------------------------------------------

      integer wrkr_sleep_event
      integer wrkr_run_aces2_event 
      parameter (wrkr_sleep_event     = 4445)
      parameter (wrkr_run_aces2_event = 4446)

c-------------------------------------------------------------------------
c   Message lengths
c-------------------------------------------------------------------------

      integer len_wrkr_integral_request
      parameter (len_wrkr_integral_request = 11)
      integer len_wrkr_contraction_request
      parameter (len_wrkr_contraction_request = 14)
