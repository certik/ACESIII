
/*
 * This routine sets n bytes of string b to the lowest byte of f_int c.
 */

#include <string.h>

#include "f_types.h"

void c_memset (void * b, f_int * c, f_int * n)
{ memset(b,(int)*c,(size_t)*n); return; }

void c_memset_(void * b, f_int * c, f_int * n)
{ memset(b,(int)*c,(size_t)*n); return; }

void C_MEMSET (void * b, f_int * c, f_int * n)
{ memset(b,(int)*c,(size_t)*n); return; }

