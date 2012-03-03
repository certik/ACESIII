
/*
 * This routine returns the elapsed real time in whole seconds and microseconds
 * since 1 January 1970 GMT.
 */

#include <stdio.h>
#include <sys/time.h>

#include "f77_name.h"
#include "f_types.h"

void
F77_NAME(c_gtod,C_GTOD)
(f_int * rtime_sec, f_int * rtime_usec)
{
    struct timeval  tp;
    struct timezone junk;
    if (gettimeofday(&tp,&junk))
        printf("@C_GTOD: gettimeofday failed -> "
               "wallclock timing stats are meaningless\n");
    else
    {
        *rtime_sec  = tp.tv_sec;  /* seconds */
        *rtime_usec = tp.tv_usec; /* microseconds */
    }
    return;
}

