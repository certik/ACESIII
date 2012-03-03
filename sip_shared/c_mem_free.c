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
#include <stdlib.h>
#include "f77_name.h"
#include "f_types.h"

/******************************************************************************
 *
 *   FORTRAN-callable memory free routine.
 *
 ****************************************************************************/
void F77_NAME(c_mem_free, C_MEM_FREE)(void *p)
{
   free(p);
}
