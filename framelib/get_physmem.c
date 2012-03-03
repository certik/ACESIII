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
#include "f77_name.h"
#include "f_types.h"

#ifdef AIX
#include <sys/systemcfg.h>    
#endif

#ifdef SGI
#include <sys/types.h>
#include <sys/sysmp.h>
#endif

#ifdef HP
#include <mach/mach_traps.h>
#include <mach/vm_statistics.h>
#endif

/**************************************************************************
 *
 *   Returns the amount of physical memory in megabytes.
 *
 **************************************************************************/

f_int F77_NAME(get_physmem, GET_PHYSMEM) ()
{
   f_int megabytes;
   size_t malloc_size;
   void *ptr;

#ifdef AIX
   extern struct _system_configuration;
   megabytes = (f_int) (_system_configuration.physmem >> 20); 
   if (megabytes > 2048) megabytes = 2048;
#endif

#ifdef SGI

#define pagetok(pages) ((((uint64_t) pages) * pagesize) >> 20)

/*********************************************************************
 *  Code to get the available memory on a system running IRIX6.2 
 *  or above.
 *********************************************************************/

   struct rminfo   realmem;
   int pagesize;

   if (sysmp(MP_SAGET, MPSA_RMINFO, &realmem, sizeof(realmem)) == -1) {
      perror("sysmp(MP_SAGET,MPSA_RMINFO, ...)");
      F77_NAME(c_exit, C_EXIT)(1);
   }

   pagesize  = getpagesize();
   megabytes = pagetok(realmem.physmem);
#endif

#ifdef HP
/************************************************************************
 *  We get the total amount of memory from the virtual memory statistics
 *  (ala the "top" command).  The total amount of memory is the sum of
 *  the free_count, active_count, inactive_count, and wire_count fields.
 *  This total must be converted from pages into megabytes.
 ************************************************************************/
/*    ###### Does not work on emerald ############
 
   long physmem;

   vm_statistics_data_t vmstats;
   int pagesize;

   pagesize = getpagesize();
   (void) vm_statistics(task_self(),&vmstats);
   physmem = pagesize * (vmstats.free_count + vmstats.active_count +
             vmstats.inactive_count + vmstats.wire_count);
    ######################################### */
   
   /* megabytes = physmem >> 20; */  /* convert to Mbytes. */
   megabytes = 1024;
#endif

#ifdef __crayx1
   megabytes = 800;   /* Hard-coded to 1 gbyte per processor on X1 */

   megabytes = F77_NAME(check_malloc, CHECK_MALLOC)();
   return (megabytes);
#endif

/* Find the maximum memory that we can actually perform a succesful malloc */

   do 
   {
      malloc_size = megabytes << 20;
      ptr = malloc(malloc_size);

      megabytes--;
   } while ( (ptr == NULL) && (megabytes > 0));

   free(ptr);
   megabytes = malloc_size >> 20;

   /* The value of "megabytes" is the most we can actually malloc.  
      Return 95% of this value to allow the operating system room for
      paging. */

   megabytes = 0.95* megabytes;
   return (megabytes);
}
