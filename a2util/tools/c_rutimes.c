
#if !defined(_CRAY_PVP) && !defined(_CRAY_MPP) && defined(__ansi_c)

/*
 * This routine returns the elapsed user and system times in whole seconds and
 * microseconds since starting the process. The infamous "CPU" time is the sum
 * of user and system times.
 */

/*
 * PORTING:
 *  - older Cray platforms (viz. SV1 generation) do not have getrusage, so
 *    rename this file to c_rutimes.F
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#define RUSAGE_SELF 0
#define RUSAGE_CHILDREN -1

#include "f77_name.h"
#include "f_types.h"

void
F77_NAME(c_rutimes,C_RUTIMES)
(f_int * utime_sec, f_int * utime_usec,
 f_int * stime_sec, f_int * stime_usec)
{
    struct rusage rc;
    if (getrusage(RUSAGE_SELF,&rc))
        printf("@C_RUTIMES: getrusage failed -> "
               "CPU timing stats are meaningless\n");
    else
    {
        *utime_sec  = rc.ru_utime.tv_sec;
        *utime_usec = rc.ru_utime.tv_usec;
        *stime_sec  = rc.ru_stime.tv_sec;
        *stime_usec = rc.ru_stime.tv_usec;
    }
    return;
}

#else /* !_CRAY_NV1 && __fortran77 */

      subroutine c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      implicit none

      integer utime_sec, utime_usec, stime_sec, stime_usec
      double precision utime

      call tsecnd(utime)
      utime_sec  = utime
      utime_usec = (utime-utime_sec)*1.d6
      stime_sec  = 0
      stime_usec = 0

      return
c     end subroutine c_rutimes
      end

#endif

