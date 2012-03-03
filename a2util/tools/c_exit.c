
/*
 * This routine exits with status.
 */

#include <stdlib.h>

#include "f_types.h"

void c_exit (f_int * status)
{ exit((int)*status); }

void c_exit_(f_int * status)
{ exit((int)*status); }

void C_EXIT (f_int * status)
{ exit((int)*status); }

