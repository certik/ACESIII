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

void F77_NAME(pack_pardo_timer, PACK_PARDO_TIMER) (f_int *op, 
                                          f_int *key1, f_int *key2)
{
   f_int x;

   x = (*key1 << 16) | (*key2);
   *op = x;
}

void F77_NAME(unpack_pardo_timer, UNPACK_PARDO_TIMER) (f_int *op, 
                                         f_int *key1, f_int *key2)
{
   f_int x;
   f_int mask;

   x = *op;
   *key1 = (x >>16);
   *key2 = (*key1 << 16) ^ x;
}
