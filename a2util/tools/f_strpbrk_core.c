
/*
 * This routine returns the FORTRAN INDEX in sz of any single char
 * appearing in charset. Both sz and charset must be NULL-terminated.
 */

#include <string.h>

#include "f77_name.h"
#include "f_types.h"

f_int
F77_NAME(f_strpbrk_core,F_STRPBRK_CORE)
(const char * sz, const char * charset)
{
    char * cz = strpbrk(sz,charset);
    if (cz) return (f_int)(1+cz-sz);
    else    return (f_int)0;
}

