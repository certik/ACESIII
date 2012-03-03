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
#include "f77_name.h"
#include "f_types.h"

void F77_NAME(convert_to_base_index, CONVERT_TO_BASE_INDEX) 
(char *anchor, long long *offset, void *base, f_int *type, long long *index)
{
   f_int *iptr;
   double *dptr;
   float *fptr;

   if (*type == 1) 
   {
      iptr = (f_int *) &anchor[*offset];
      *index = 1 + iptr - (f_int *) base;
      return;
   }

   if (*type == 2) 
   {
      fptr = (float *) &anchor[*offset];
      *index = 1 + fptr - (float *) base;
      return;
   }

   if (*type == 3) 
   {
      dptr = (double *) &anchor[*offset];
      *index = 1 + dptr - (double *) base;
      return;
   }
}

