
c This routine is the common exit handler for ACES.

      subroutine aces_exit(iStat)
      implicit none
      integer iStat, iErr
#include "mpif.h"
#include "parallel_aces.com"
      if (irank.eq.0) then
         if (iStat.eq.0) then
            print *,
     &          'The ACES Member Executable has completed successfully.'
         else
            print *, '@ACES_EXIT: An ACES error has occurred.'
            print *, '            return status = ',iStat
         end if
      end if
      if (nprocs.eq.1) then
#ifdef _UNICOS /* CRAY intrinsics */
         if (iStat.ne.0) then
#ifndef _CRAY_NV1
            call tracebk
            call trbk
#endif
            call abort
         end if
#endif /* _UNICOS */
         call c_exit(iStat)
      else
#ifdef _HAVE_MPI
         call MPW_Abort(MPI_COMM_WORLD,iStat,iErr)
#else
c      o this branch makes no sense (nprocs!=1 w/o MPI), but it hides
c        mpw from parallel builds that use serial ACES libraries
         call c_exit(iStat)
#endif
      end if
      end

