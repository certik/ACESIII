
/*
 * This routine mimics (blas)dcopy for integers.
 */

#include <string.h>
#include <stdlib.h> /* for labs */

#include "f77_name.h"
#include "f_types.h"

#define _ABORT_ON_OVERLAP

void
F77_NAME(icopy,ICOPY)
(f_int * len, f_int * src, f_int * inc_src,
              f_int * dst, f_int * inc_dst)
{
    long i;

    if ((*len <= 0) || !*inc_src || !*inc_dst) return;

    if ( (*inc_src==1) && (*inc_dst==1) )
    {
        if (dst==src) return;
#ifdef _DEBUG
     /* complain about overlapping data */
        i = labs(src-dst);
        if (i<*len)
        {
            printf("\n"
                   "@ICOPY: WARNING - the source and destination data overlap\n"
                   "        src address: %p (%li)\n"
                   "        dst address: %p (%li)\n"
                   "        difference:  %li f_ints\n"
                   "        data length: %li f_ints\n\n",
                   src,(long)src,
                   dst,(long)dst,
                   i,
                   (long)*len
                  );
#ifdef _ABORT_ON_OVERLAP
            abort();
        }
#endif
#endif
        memmove(dst,src,(size_t)(*len*sizeof(f_int)));
    }
    else
    {
     /* YAU : old
        for (i=0; i<*len; i++) *(dst+i*(*inc_dst)) = *(src+i*(*inc_src));
        YAU : new
              This is an induction variable optimization resulting in
              two strength reductions. */
        f_int *pDst = dst, *pSrc = src;
        for (i=0; i<*len; i++)
        {
            *pDst = *pSrc;
            pDst += *inc_dst;
            pSrc += *inc_src;
        }
     /* YAU : end */
    }

    return;
}

