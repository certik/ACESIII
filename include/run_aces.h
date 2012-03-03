c-----------------------------------------------------------------------
c   Parameter common block for the "run_aces" mission.
c-----------------------------------------------------------------------

      common /run_aces/first_run_aces, aces_dir, last_run_aces
      integer first_run_aces, last_run_aces
      character*80 aces_dir
