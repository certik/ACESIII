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
 * This routine frees a block of memory allocated with malloc.
 *
 * Mixing and matching malloc/free with FORTRAN programs is tricky.
 * Do we free the address or the value at the address? If the latter,
 * then is it an f_int or an f_adr type? Since this is supposed to be
 * used in conjunction with c_malloc (or c_sbrk), it will free the f_adr
 * value at the address. Observe:
 *
 * #include "f_types.h"
 *       integer bytes
 *       F_ADR z, c_malloc
 *       z = c_malloc(bytes)
 *       call c_free(z)
 */

#include <stdio.h>
#include <stdlib.h>
#include "f77_name.h"
#include "f_types.h"

void F77_NAME(c_free, C_FREE)(f_adr * n)
{ printf("C_FREE: n = %p\n",*n); free((void *)*n); return; }

