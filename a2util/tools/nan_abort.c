
/*
 * This routine aborts if the double argument is NaN.
 */

#include <stdlib.h>

#include "f77_name.h"

void
F77_NAME(nan_abort,NAN_ABORT)
(unsigned char * p)
{ if ( (*(p+0)==0xff) && (*(p+1)==0xff) &&
       (*(p+2)==0xff) && (*(p+3)==0xff) &&
       (*(p+4)==0xff) && (*(p+5)==0xff) &&
       (*(p+6)==0xff) && (*(p+7)==0xff)
     ) abort();
  return;
}

