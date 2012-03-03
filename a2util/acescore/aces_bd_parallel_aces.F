
c This routine defines the default values for the common blocks in
c parallel_aces.com.

      blockdata aces_bd_parallel_aces
      implicit none
#define NO_EXTERNAL
#include "parallel_aces.com"
      data nprocs,irank /1,0/
      data szcpuname /'localhost'/
      data  icpuname /9/
      end

