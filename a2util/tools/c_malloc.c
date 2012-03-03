
/*
 * This routine allocates a block of memory of size n bytes.
 */

#include <stdlib.h>

#include "f_types.h"

f_adr c_malloc (f_int * n)
{ return (f_adr)malloc(*n); }

f_adr c_malloc_(f_int * n)
{ return (f_adr)malloc(*n); }

f_adr C_MALLOC (f_int * n)
{ return (f_adr)malloc(*n); }

