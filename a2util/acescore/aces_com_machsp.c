
/*
 * This routine initializes the machsp common block.
 */

#include <stdio.h>
#include <stdlib.h>

#include "f77_name.h"
#include "f_types.h"

#include "machsp.com"

void
F77_NAME(aces_com_machsp,ACES_COM_MACHSP)
()
{
    f_int iTmp;
    const f_int  i = sizeof(i);
    const double d = sizeof(d);
#ifdef _ASSERT
#ifdef F_64BIT
    if (i != 8)
    {
        printf("@ACES_COM_MACHSP: Assertion failed.\n"
               "                  sizeof(f_int) = %i\n"
               "                  F_64BIT       = 1\n",i);
        abort();
    }
#else
    if (i != 4)
    {
        printf("@ACES_COM_MACHSP: Assertion failed.\n"
               "                  sizeof(f_int) = %i\n"
               "                  F_64BIT       = 0\n",i);
        abort();
    }
#endif
#endif /* _ASSERT */
    f_machsp.iintln = sizeof(f_int);	/* bytes per integer */
    f_machsp.ifltln = sizeof(double);	/* bytes per double */
    f_machsp.iintfp = d/i;		/* integers per double */
    iTmp = i<<1;
    f_machsp.ialone = (1<<iTmp)-1;	/* bitmask for 1/4 of an integer */
    f_machsp.ibitwd = iTmp;		/* bits per 1/4 of an integer */
    return;
}

