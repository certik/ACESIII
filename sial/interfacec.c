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
#include "parser_interface.h"
/*#include "testf.h"*/
#include "f77_name.h"
#include "f_types.h"

extern int mx_array_dim;

int create_optable_c(int nentries)
{
    int rc;
    f_int f_nentries = nentries;
    rc = (int) F77_NAME(create_optable, CREATE_OPTABLE) (&f_nentries);
    return rc;
}

int add_optable_c(int opcode, 
                    int op1_array,
                    int op2_array,
                    int result_array,
                    int indarray[mx_array_dim],
                    int user_sub_index,
                    int lineno)
{
   int rc;
   int i;

   f_int findarray[mx_array_dim];
   f_int f_opcode = opcode;
   f_int f_op1_array = op1_array;
   f_int f_op2_array = op2_array;
   f_int f_result_array = result_array;
   f_int f_user_sub_index = user_sub_index;
   f_int f_lineno = lineno;

   if (opcode == 159 || opcode == 160 || opcode == 110)     /* allocate or create */
   {
      /* Allocate, deallocate, or create instruction requires raw index array, 
         do not add 1 */

      for (i = 0; i < mx_array_dim; i++)
         findarray[i] = indarray[i];
   }
   else
   {
      /* Normal case: Adjust indices for Fortran indexing */

      for (i = 0; i < mx_array_dim; i++)
         findarray[i] = indarray[i] + 1;
   }

   rc = (f_int) F77_NAME(add_optable, ADD_OPTABLE) (&f_opcode,
                                                    &f_op1_array,
                                                    &f_op2_array,
                                                    &f_result_array,
                                                    findarray,
                                                    &f_user_sub_index,
                                                    &f_lineno);

   return rc;
}

