/*
*  Copyright (c) 2003-2010 University of Florida
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  The GNU General Public License is included in this distribution
*  in the file COPYRIGHT.
*/ 

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

