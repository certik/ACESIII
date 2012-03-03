
/*
 * This routine zeroes a double array of length elements.
 */

#include "f77_name.h"
#include "f_types.h"

void
F77_NAME(dzero,DZERO)
(double * dArr, f_int * length)
{
    if (*length>0)
    {
        long i;
        for (i=0;i<*length;i++) *(dArr+i) = 0.0;
    }
    return;
}

