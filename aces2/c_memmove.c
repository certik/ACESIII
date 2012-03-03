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

