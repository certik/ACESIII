#include <stdio.h>
#include "f_types.h"

/* Declarations of blkmgr interface routines. */

double *c_get_block_addr(int arrayno, int blocknumber, int blkndx);
double *c_get_block_data_addr(int arrayno, int blocknumber, int blkndx);
int c_allocate_block(int arrayno, int blocknumber);
int c_get_block_number(int arrayno, int blocknumber, f_int *array_table,
         f_int *block_map_table);
