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
#include <ctype.h>
#include <string.h>
#include "f77_name.h"
#include "f_types.h"

void F77_NAME(c_line_parse, C_LINE_PARSE) 
(char *line, char *token, f_int *flag, f_int *f_ival, double *rval, char *cval)
{
    int i;
    int ival;

    /* Scan for the token.  If not found, return. */
 
    for (i = 0; i < strlen(line); i++)
       if (line[i] == *token) 
       {
          i++; 
          while (i < strlen(line))
          {
             if (!isspace(line[i])) break; 
             i++;
          }
          
          switch (*flag)
          {
             case 1:   sscanf(&line[i], "%d", &ival);  *f_ival = ival; /* integer */
             case 2:   sscanf(&line[i], "%lf", rval);  /* double */
             case 3:   strcpy(cval, &line[i]);         /* character */
          } 
       } 
}
           
void F77_NAME(c_decode_csv_integer, C_DECODE_CSV_INTEGER)(char *s, f_int *x, f_int *nvals)
{
   /* Fortran-callable routine to decode a comma- or blank-separated string into
      an integer array. The data is returned in x, and the number of decoded
      values is returned in nvals. */

   int i, n, nn, istart, iend;
   char cstr[100];
   
   n = 0;
   istart = 0;

   while (istart < strlen(s))
   {
      /* Skip leading whitespace. */

      while (istart < strlen(s))
      {
         if (!isspace(s[istart])) break;
         istart++;
      }

      /* Find trailing separator char (either blank or comma */

      iend = istart;
      for (i = iend+1; i < strlen(s); i++)
      {
         if (s[i] == ' ' || s[i] == ',') break;
         iend++; 
      }

      /* Decode the next value */
    
      strncpy(cstr, &s[istart], iend-istart+1);
      nn = sscanf(cstr, "%d", &x[n]); 
      if (nn > 0) 
      {
         n++;

         istart = i + 1;
      }
      else
      {
         break;
      }
      
   }  /* while (istart < strlen(s) */
   *nvals = n;
}
 

