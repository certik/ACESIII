
/*
 * This routine moves/copies len bytes from src to dst.
 */

#include <string.h>

#include "f_types.h"

void c_memmove (void * dst, void * src, f_int * len)
{ memmove(dst,src,(size_t)*len); return; }

void c_memmove_(void * dst, void * src, f_int * len)
{ memmove(dst,src,(size_t)*len); return; }

void C_MEMMOVE (void * dst, void * src, f_int * len)
{ memmove(dst,src,(size_t)*len); return; }

