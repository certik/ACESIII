#ifndef _PARALLEL_ACES_COM_
#define _PARALLEL_ACES_COM_
c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

#ifndef NO_EXTERNAL
      external aces_bd_parallel_aces
#endif

#include "parallel_aces.h"

      integer                nprocs, irank, icpuname

      character*(_MPI_MAX_PROCESSOR_NAME) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end
#endif /* _PARALLEL_ACES_COM_ */
