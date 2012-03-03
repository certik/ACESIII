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

