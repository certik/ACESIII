
/*
 * This routine returns the exit value of the system call on *str.
 */

#include <stdlib.h>

#include "f_types.h"

f_int c_system (const char * str)
{ return (f_int)system(str); }

f_int c_system_(const char * str)
{ return (f_int)system(str); }

f_int C_SYSTEM (const char * str)
{ return (f_int)system(str); }

