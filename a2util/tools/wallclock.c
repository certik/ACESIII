
/*
 * This routine returns the current system time since the Epoch
 * (00:00:00 UTC, 1 January 1970).
 */

#include <time.h>

#include "f77_name.h"
#include "f_types.h"

void
F77_NAME(wallclock,WALLCLOCK)
(f_int * year, f_int * mon, f_int * mday,
 f_int * hour, f_int * min, f_int * sec)
{
    struct tm  now;
    struct tm *junk;

    time_t t = time(0);
    junk = localtime_r(&t,&now);

    if (now.tm_year < 1900)
        *year = now.tm_year + 1900;
    else
        *year = now.tm_year;
    *mon  = now.tm_mon + 1;
    *mday = now.tm_mday;
    *hour = now.tm_hour;
    *min  = now.tm_min;
    *sec  = now.tm_sec;

    return;
}

