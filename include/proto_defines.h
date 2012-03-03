c------------------------------------------------------------------------
c   Defined constants for use in the proto program.
c------------------------------------------------------------------------

       integer master_status
       integer worker_status
       integer manager_status
       integer both_status
       integer master_worker_status
       integer unassigned_status
       integer terminate_state
       integer wrkr_terminate_state
       integer mgr_terminate_state
       integer ready_state
       integer mgr_ready_state
       integer wrkr_ready_state
       integer uninitialized_state
       integer mstr_error_state

       parameter (worker_status  = 1)
       parameter (manager_status = 2)
       parameter (both_status    = 3)
       parameter (master_status  = 4)
       parameter (master_worker_status = 5)
       parameter (unassigned_status = 0)

       parameter (terminate_state = 101)   ! indicates worker termination
       parameter (ready_state     = 102)   ! all ready to handle events
       parameter (wrkr_terminate_state = 103) ! worker has terminated
       parameter (mgr_ready_state = 104)   ! mgr ready to handle events
       parameter (wrkr_ready_state = 105)   ! wrkr ready to handle events
       parameter (mgr_terminate_state = 106) ! manager has terminated
       parameter (uninitialized_state = 203)  
       parameter (mstr_error_state = 204)  
