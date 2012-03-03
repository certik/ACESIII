c--------------------------------------------------------------------------
c   Tracelevel flags
c--------------------------------------------------------------------------

      integer instruction_trace
      parameter (instruction_trace = 1)
      integer contraction_trace
      parameter (contraction_trace = 2)
      integer proc_trace
      parameter (proc_trace = 4)
   
      common /trace/tracelevel, current_op, current_line,trace,
     *              call_marker, dryrun, simulator, pardo_timer,
     *              pardo_block_wait_timer
      integer tracelevel
      integer current_op
      integer current_line
      integer call_marker
      integer pardo_timer, pardo_block_wait_timer
      logical trace, dryrun, simulator
