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
#include "f_types.h"
#include "f77_name.h"

#define MAX_USER_SUB 500
#define MAX_USER_SUB_NAME 33

int nuser_sub = 0;
void *user_sub[MAX_USER_SUB];
char user_sub_names[MAX_USER_SUB][MAX_USER_SUB_NAME];

extern void F77_NAME(abort_job,ABORT_JOB)();

void F77_NAME(clear_user_sub, CLEAR_USER_SUB)()
{
     nuser_sub = 0;
}

f_int F77_NAME(load_user_sub, LOAD_USER_SUB)(char *sub_name,  void (*fp)() )
{
   /* Stores the function pointer passed in as argument in the user_sub table,
      returning a handle to use when executing the function. */

   int handle, ierr;

   nuser_sub++;
   if (nuser_sub > MAX_USER_SUB) 
   {
      printf("user_sub table is full.  The current limit is %d.\n", 
             MAX_USER_SUB);
      ierr = 1;
      F77_NAME(abort_job, ABORT_JOB)();
   }
 
   handle = nuser_sub;
   user_sub[nuser_sub-1] = (void *)fp;
   if (strlen(sub_name) <= MAX_USER_SUB_NAME)
      strcpy(user_sub_names[nuser_sub-1], sub_name);
   else
   {
      printf("User subroutine name %s greater than %d chars.\n", 
             sub_name, MAX_USER_SUB_NAME);
      ierr = 2;
      F77_NAME(abort_job, ABORT_JOB)();
   }

   return ( (f_int) handle);
}

void F77_NAME(exec_user_sub, EXEC_USER_SUB) (f_int *handle, 
      f_int *array_table, f_int *narray_table, 
      f_int *index_table, f_int *nindex_table, 
      f_int *segment_table, f_int *nsegment_table,
      f_int *block_map_table, f_int *nblock_map_table,
      double *scalar_table, f_int *nscalar_table,
      long long *address_table,
      f_int *op ) 
{
   /* Executes the subroutine stored in the user_sub table, referenced by
      handle. */

   int ierr;
   void (*fp)(f_int *, f_int *, f_int *, f_int *, f_int *, f_int *,
              f_int *, f_int *, double *, f_int *, long long *, f_int *);

   if (*handle > 0 && *handle <= MAX_USER_SUB) 
   {
      fp = (void(*)())user_sub[*handle-1];
      fp (array_table, narray_table, index_table, nindex_table, 
          segment_table, nsegment_table, block_map_table,
          nblock_map_table, scalar_table, nscalar_table,
          address_table, op);
   }
   else
   {
      printf("Invalid handle in exec_user_sub: %d.\n", *handle);
      ierr = 2;
      F77_NAME(abort_job, ABORT_JOB)();
   }
}

int c_get_subroutine_handle(char *sub_name)
{
   /* Returns the key associated with a pre-defined subroutine name. */
   /* If the subroutine is not found, it returns -1. */

   int i;

   for (i = 0; i < nuser_sub; i++)
   {
      if (!strcmp(user_sub_names[i], sub_name) )
         return (i+1);
   }

   return -1;
}

f_int F77_NAME(get_subroutine_handle, GET_SUBROUTINE_HANDLE) (char *sub_name)
{
   return c_get_subroutine_handle(sub_name);
}
