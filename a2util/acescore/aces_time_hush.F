
c This routine switches off the timing statistics that are printed by
c aces_fin.

      subroutine aces_time_hush
      implicit none
#include "aces_time.com"
      ame_timed = .false.
      return
      end

