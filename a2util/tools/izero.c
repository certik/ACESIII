
/*
 * This routine zeroes a f_int array of length elements.
 */

#include "f77_name.h"
#include "f_types.h"

void
F77_NAME(izero,IZERO)
(f_int * lArr, f_int * length)
{
    if (*length>0)
    {
        long i;
        for (i=0;i<*length;i++) *(lArr+i) = 0;
    }
    return;
}

