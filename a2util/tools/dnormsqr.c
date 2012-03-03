
/*
 * This routine returns the dot product of a vector with itself, i.e.
 * the square of the norm. (This routine is faster than ddot in principle
 * because it knows the two arrays are aliased.)
 */

#include "f77_name.h"
#include "f_types.h"

double
F77_NAME(dnormsqr,DNORMSQR)
(f_int * len, double * dVec, f_int * incr)
{
    double dTmp = 0.0;
    if ((*len > 0) && *incr)
    {
        if (*incr==1)
        {
            long i;
            for (i=0; i<*len; i++)
                dTmp += (*(dVec+i))*(*(dVec+i));
        }
        else
        {
            double *pd = dVec;
            long i;
            for (i=0; i<*len; i++)
            {
                dTmp += (*pd)*(*pd);
                pd   += *incr;
            }
        }
    }
    return dTmp;
}

