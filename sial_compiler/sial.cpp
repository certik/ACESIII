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
#include "interpret.h"
#include "parser_interface.h"
#include "f77_name.h"
#include "f_types.h"
#include "string.h"

extern "C" {
void F77_NAME(load_pre_defined_routines_compile_time, 
              LOAD_PRE_DEFINED_ROUTINES_COMPILE_TIME) ();

void F77_NAME(mem_alloc_init, MEM_ALLOC_INIT)(f_int *mem, f_int *sheap_flag,
                                              f_int *ierr);

f_int F77_NAME(get_max_array_dimension, GET_MAX_ARRAY_DIMENSION)();

void F77_NAME(abort_job, ABORT_JOB)()
{
   exit(1);
}

#ifdef XT3
void F77_NAME(pghpf_0l, PGHPF_0L)()
{}
#endif

};

int mx_array_dim;

int main (int argc, const char* const* argv)
{
    int argnumber=1;
    int printindextable=0;
    int printarraytable=0;
    int printinstructiontable=0;
    int printparse=0;

    int hasinfilename=0;
    char infilename[240];
    int length=0;
    int len;

    int defaultoutfilename=1;
    char outfilename[240];
    f_int memory_required;
    f_int sheap_flag = 0;
    f_int ierr;
    f_int f_mx_array_dim;

    if (argc==1)
    {
        printf("No input file.\n");
        printf("Please use the format:\n");
        printf("  %s sialfilename.sial <-o objfilename> <-i> <-a> <-s>\n",
            argv[0]);
        /*printf("-l: lexical analysis results\n");*/
        printf("-i: print out index table results\n");
        printf("-a: print out array table results\n");
        printf("-s: print out instruction table results\n");
        return 1;
    }

    /* Set up pre-defined subroutine names for later syntax checking */

    F77_NAME(load_pre_defined_routines_compile_time, 
             LOAD_PRE_DEFINED_ROUTINES_COMPILE_TIME) ();

    /* Get the max. dimensions supported by the Fortran-based code. */

    f_mx_array_dim = F77_NAME(get_max_array_dimension, 
                              GET_MAX_ARRAY_DIMENSION)();
    mx_array_dim = f_mx_array_dim;
 
    /* Initialize memory subsystem used by Fortran routines. */

    memory_required = 300;
    F77_NAME(mem_alloc_init, MEM_ALLOC_INIT)(&memory_required,
                                             &sheap_flag, &ierr);
    while (argnumber<argc)
    {
        if (!strncmp(argv[argnumber], "-o", 2)&&argnumber+1<argc)
        {
            argnumber++;
            defaultoutfilename=0;
            len=strlen(argv[argnumber]);
            strncpy(outfilename, argv[argnumber], len);
            outfilename[len]='\0';
            argnumber++;
        }
        else if (!strncmp(argv[argnumber], "-i", 2))
        {
            argnumber++;
            printindextable=1;
        }
        else if (!strncmp(argv[argnumber], "-a", 2))
        {
            argnumber++;
            printarraytable=1;
        }
        else if (!strncmp(argv[argnumber], "-s", 2))
        {
            argnumber++;
            printinstructiontable=1;
        }
        else if (!strncmp(argv[argnumber], "-l", 2))
        {
            argnumber++;
            printparse=1;
        }
        else
        {
            length=strlen(argv[argnumber]);
            if (length>5&&!strncmp(".sial", (const char*)
                                     (&(argv[argnumber][length-5])), 5))
            {
                strncpy(infilename, argv[argnumber], length);
                infilename[length]='\0';
                hasinfilename=1;
            }
            else
            {
                printf("unknown argument: %s \n", argv[argnumber]);
                return 1;
            }
            argnumber++;
        }
    }/* end of while*/

    if (!hasinfilename)
    {
        printf("No input file.\n");
        return 1;
    }

    if (defaultoutfilename)
    {
        strncpy(outfilename, infilename, length-5);
        outfilename[length-5]='\0';
        strncat(outfilename, ".sio", 4);
        outfilename[length-5+4]='\0';
    }

    printf("input  file: %s\noutput file: %s\n ", 
        infilename, outfilename);

    if (interpret(infilename, printparse,printindextable, printarraytable,
                    printinstructiontable)==0)
    {
        printf("Writing tables to %s\n",outfilename);
        F77_NAME(write_tables, WRITE_TABLES)(outfilename);
        printf("\nSIO (Super Instruction Object) file is saved to \"%s\".\n",
                 outfilename);
    }
    return 0;
}
