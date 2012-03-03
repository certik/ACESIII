
/*
 * This routine removes a file or directory.
 */

#include <stdio.h>

#include "f_types.h"

f_int f_remove_core (const char * file)
{ return (f_int)remove(file); }

f_int f_remove_core_(const char * file)
{ return (f_int)remove(file); }

f_int F_REMOVE_CORE (const char * file)
{ return (f_int)remove(file); }

