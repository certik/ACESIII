
/*
 * This routine frees a block of memory allocated with malloc.
 *
 * Mixing and matching malloc/free with FORTRAN programs is tricky.
 * Do we free the address or the value at the address? If the latter,
 * then is it an f_int or an f_adr type? Since this is supposed to be
 * used in conjunction with c_malloc (or c_sbrk), it will free the f_adr
 * value at the address. Observe:
 *
 * #include "f_types.h"
 *       integer bytes
 *       F_ADR z, c_malloc
 *       z = c_malloc(bytes)
 *       call c_free(z)
 */

#include <stdlib.h>

#include "f_types.h"

void c_free (f_adr * n)
{ free((void *)*n); return; }

void c_free_(f_adr * n)
{ free((void *)*n); return; }

void C_FREE (f_adr * n)
{ free((void *)*n); return; }

