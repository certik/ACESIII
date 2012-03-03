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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "f77_name.h"
#include "f_types.h"

void F77_NAME(generate_scratch_filename, 
              GENERATE_SCRATCH_FILENAME) (char *fn, char *extension)
{
   /*   Generates a unique filename using the form SCRxxxxxx.extension.
        The filename is returned in the "fn" argument.                  */
   char template[] = "SCR_XXXXXX";
   char *cp;

   cp = mktemp(template);
   strcpy(fn, cp);
   strcat(fn,".");
   strcat(fn,extension);
}
