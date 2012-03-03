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
#include <stdlib.h>
#include "f_types.h"
#include "f77_name.h"
#ifdef HP
#include <stddef.h>
#endif

long long F77_NAME(unmask_addr, UNMASK_ADDR)(long long *addr)
{
    /* Masks off extraneous bits used on addresses on Altix and 
       other systems.  The unmasked address is returned by the
       function.  This allows Fortran address arithmetic to be 
       performed on these systems. */

   long long raddr, mask;

#ifdef ALTIX
    /* mask = 0x00ffffffffffffff; */
    mask = 0x0000000fffffffff;
    raddr = mask & (*addr);
#else
   /* No masking.  Pass back the raw address. */

   raddr = *addr;
#endif

   return raddr;
}


