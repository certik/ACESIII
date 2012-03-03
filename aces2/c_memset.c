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

