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

