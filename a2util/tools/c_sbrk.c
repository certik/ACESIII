
/*
 * This routine increments the program break point by n bytes.
 * YAU: I'm hoping that C compilers will automatically recast the "*n"
 *      argument since different OSes take different types.
 */

#include <unistd.h>

#include "f_types.h"

f_adr c_sbrk (f_int * n)
{ return (f_adr)sbrk(*n); }

f_adr c_sbrk_(f_int * n)
{ return (f_adr)sbrk(*n); }

f_adr C_SBRK (f_int * n)
{ return (f_adr)sbrk(*n); }

