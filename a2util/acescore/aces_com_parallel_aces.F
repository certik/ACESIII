
c This routine initializes the parallel_aces common block.

      subroutine aces_com_parallel_aces
      implicit none

#include "mpif.h"

c INTERNAL VARIABLES
      integer iErr
      logical bMPIUp

c COMMON BLOCKS
#include "parallel_aces.com"

c ----------------------------------------------------------------------

#ifdef _HAVE_MPI
c   o return if MPI is not initialized
      call MPW_Initialized(bMPIUp,iErr)
      if (iErr.ne.0) then
         print *, '@ACES_COM_PARALLEL_ACES: ',
     &            'unable to determine MPI initialization status'
         call aces_exit(1)
      end if
      if (.not.bMPIUp) return

c   o update the number of MPI processes
      call MPW_Comm_size(MPI_COMM_WORLD,nprocs,iErr)
      if (iErr.ne.0) then
         print *, '@ACES_COM_PARALLEL_ACES: ',
     &            'unable to determine MPI comm size'
         call aces_exit(1)
      end if

c   o update the process rank
      call MPW_Comm_rank(MPI_COMM_WORLD,irank,iErr)
      if (iErr.ne.0) then
         print *, '@ACES_COM_PARALLEL_ACES: ',
     &            'unable to determine MPI comm rank'
         call aces_exit(1)
      end if

c   o bomb if the hard-coded limit has been exceeded
      if (nprocs.gt._PA_MAX_PROCS) then
         if (irank.eq.0) then
         print *, '@ACES_COM_PARALLEL_ACES: The maximum number of ',
     &            'processes allowed by the'
         print *, '                         ACES subsystem has been ',
     &            'exceeded. Either reduce'
         print *, '                         the requested number to at',
     &            ' most ',_PA_MAX_PROCS
         print *, '                         or recompile the source ',
     &            'with a higher limit.'
         end if
c      o make sure the root has a chance to print the error
         call MPW_Barrier(MPI_COMM_WORLD,iErr)
         call aces_exit(1)
      end if

c   o update the machine name (essentially `uname -n`)
      call MPW_Get_processor_name(szcpuname,icpuname,iErr)
      if (iErr.ne.0) then
         print *, '@ACES_COM_PARALLEL_ACES: ',
     &            'unable to determine MPI processor name'
         call aces_exit(1)
      end if
#endif /* _HAVE_MPI */

#ifdef _DEBUG
      if (nprocs.ne.1) then
         print *, '@ACES_COM_PARALLEL_ACES: ',
     &            'irank = ',irank,'; nprocs = ',nprocs,
     &            '; processor name = "',szcpuname(1:icpuname),'"'
      end if
#endif

      return
c     end subroutine aces_com_parallel_aces
      end

