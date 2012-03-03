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
 * This routine compares two strings up to len bytes ignoring case.
 */

#include <strings.h>

#include "f_types.h"

f_int c_strncasecmp (const char * s1, const char * s2, f_int * len)
{ return (f_int)strncasecmp(s1,s2,(size_t)*len); }

f_int c_strncasecmp_(const char * s1, const char * s2, f_int * len)
{ return (f_int)strncasecmp(s1,s2,(size_t)*len); }

f_int C_STRNCASECMP (const char * s1, const char * s2, f_int * len)
{ return (f_int)strncasecmp(s1,s2,(size_t)*len); }

