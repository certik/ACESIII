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
#ifdef HP
#include <stddef.h>
#endif

#ifdef ALTIX
#include <mpp/shmem.h>
#endif

#define MAX_MEM_ALLOC_PTRS 1000
f_int *malloced_pointers[MAX_MEM_ALLOC_PTRS];
int malloced_len[MAX_MEM_ALLOC_PTRS];
int n_malloced_pointers = 0;

void F77_NAME(malloc_wrapper, MALLOC_WRAPPER)(f_int *nwords,
                f_int *element_size, f_int *sheap_flag,
                void *x, long long *ixmem, 
                f_int *ierr)
{
   size_t nbytes = (*nwords) * (*element_size);
   size_t nb_shmalloc;
   f_int *fp, *fp2;
   long long llvar;
#ifdef HP
   ptrdiff_t offset;
#else
   long long  offset;
#endif

   if (*sheap_flag) 
   {
#ifdef USE_SHMEM
      nb_shmalloc = (size_t) nbytes;
      printf("Call shmalloc with %ld\n", nb_shmalloc);
      fp = (f_int *) shmalloc(nb_shmalloc);
      printf("Pointer from shmalloc: %x\n", fp);
#else
      printf("Error: shmalloc not implemented for this machine.\n");
      *ierr = -1; 
      return;
#endif
   }
   else
   {
      nbytes = (*nwords);
      nbytes *= (*element_size);
      nbytes += 8;

      fp = (f_int *) malloc(nbytes);
   }

   if (!fp)
   {
      printf("(SH)MALLOC ERROR: nwords %d, element size %d\n",*nwords,
              *element_size);
      *ierr = -1;
      return;
   }
  
   llvar = (long long) fp;
   if (llvar/8*8 == llvar) 
      fp2 = fp;
   else
      fp2 = fp + 1;
   
   /* Enter the memory pointer into the set of malloc'd pointers. */
   if (n_malloced_pointers >= MAX_MEM_ALLOC_PTRS)
   {
      printf("Error: Limit of %d mem_alloc pointers has been exceeded.\n",
             MAX_MEM_ALLOC_PTRS);
      *ierr = -1;
      return;
   }

   malloced_pointers[n_malloced_pointers] = fp;
   malloced_len[n_malloced_pointers] = nbytes;
   n_malloced_pointers++; 
 
    offset = ( fp2 - (f_int *) x);
   *ixmem = (long long) (( fp2 - (f_int *) x) / (*element_size/sizeof(f_int)) + 1);
    /* printf("malloc_wrapper: fp = %x %ld, address of x = %x %ld, ixmem = %ld\n",fp,fp,x,x,*ixmem);  */ 

   *ierr = 0;
}

void F77_NAME(free_mem_alloc, FREE_MEM_ALLOC) ()
{
   /* Frees all memory allocated by calls to mem_alloc */
   int i;

   for (i = 0; i < n_malloced_pointers; i++)
   {
      free( malloced_pointers[i] );
   }

   n_malloced_pointers = 0;
}

long long F77_NAME(c_loc64, C_LOC64)(char *base, long long *ix, f_int *size)
{
   long long addr;

   addr = (long long) base + (*ix-1)*(*size);
   /* printf("C_LOC64: base address %x %ld, ix %ld, size %d, addr %x %ld\n",
        base,base,*ix,*size,addr,addr); */
   return addr;
}

long long F77_NAME(c_loc64p, C_LOC64P)(char *base, long long *ix, f_int *size)
{
   long long addr;

   addr = (long long) base + (*ix-1)*(*size);
   return addr;
}


