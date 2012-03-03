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
#include <mpi.h>

/************************************************************************
 *
 *   c_blkmgr.c 
 *
 *   C interface for the blkmgr.f code.
 *
 ************************************************************************/

long long F77_NAME(get_block_index, GET_BLOCK_INDEX) (f_int *array, f_int *block,
                   f_int *stack, f_int *blkndx, double *x, f_int *data_flag);

long long F77_NAME(get_block_data_index, GET_BLOCK_DATA_INDEX) (f_int *array, 
                   f_int *block, f_int *stack, f_int *blkndx, double *x);

f_int F77_NAME(find_current_block, FIND_CURRENT_BLOCK) (f_int *x, f_int *array_table, f_int *index_table, f_int *nindex_table, f_int *segment_table, f_int *nsegment_table, f_int *blkndx);

f_int F77_NAME(get_block_number_distributed, GET_BLOCK_NUMBER_DISTRIBUTED) (f_int *array, f_int *block, f_int *array_table, f_int *block_table);

void F77_NAME(dump_block_ids, DUMP_BLOCK_IDS)();

#ifdef ALTIX
   double *dshptr;

   void F77_NAME(set_memory_pointer, SET_MEMORY_POINTER)(double *dp)
   {
      dshptr = dp;   /* Sets the memory reference pointer for Altix. */
   }
#endif

/***********************************************************************
 *   C wrapper routines that call blkmgr routines defined above.       *
 ***********************************************************************/

double *c_get_block_addr(int array, int block, int blkndx)
{
   long long ix;
   f_int farray, fblock, fblkndx, fstack, fdata_flag;
#ifdef ALTIX
   double *a = dshptr;
#else
   double a[1];
#endif
   
   farray = array;
   fblock = block;
   fblkndx = blkndx + 1; 
   fblkndx = 0;
   fdata_flag = 0;   /* Include block ID section of block. */
   ix = F77_NAME(get_block_index, GET_BLOCK_INDEX) (&farray, &fblock,
                   &fstack, &fblkndx, &(a[0]), &fdata_flag);

   ix--;      /* Compensate for FORTRAN indexing. */

   /* Return the address of a[ix] */
   return ( &(a[ix]) );
}

double *c_get_block_data_addr(int array, int block, int blkndx)
{
   long long ix;
   f_int farray, fblock, fblkndx, fstack;
#ifdef ALTIX
   double *a = dshptr;
#else
   double a[1];
#endif
   
   if (blkndx == -1) 
   {
      printf("SERVER CALLED WITH BLKNDX -1\n");
      printf("array %d array_block %d\n",array,block);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   farray = array;
   fblock = block;
   fblkndx = blkndx + 1; 
   ix = F77_NAME(get_block_data_index, GET_BLOCK_DATA_INDEX) (&farray,
                   &fblock, &fstack, &fblkndx, &(a[0]));
   
   ix--;   /* compensate for Fortran indexing */

   return ( &(a[ix]) );
}

int c_get_block_number(int array, int blocknumber, f_int *array_table, 
                       f_int *block_map_table)
{
   int rc;
   f_int arrayno, blocknum;

   arrayno = array;
   blocknum = blocknumber;
   rc = (int) F77_NAME(get_block_number_distributed, 
                       GET_BLOCK_NUMBER_DISTRIBUTED) (&arrayno, &blocknum,
                                        array_table, block_map_table);

   rc--;   /* compensate for FORTRAN block indexing */
   return (rc);
}

