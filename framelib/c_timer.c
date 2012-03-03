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
#include <time.h>
#include <sys/times.h>

void (F77_NAME(register_timer, REGISTER_TIMER)) (char *desc, f_int *type, f_int *fkey);
void (F77_NAME(update_timer, UPDATE_TIMER)) (f_int *fkey);
void (F77_NAME(timer_start, TIMER_START)) (f_int *key);
void (F77_NAME(stop_timers, STOP_TIMERS)) ();
void (F77_NAME(print_timer, PRINT_TIMER)) (char *str, f_int *key);

/*----------------------------------------------------------------------
   C interface for the timer subsystem.  Allows software written in C
   to register and update timers with the Fortran-based timer subsystem.
-------------------------------------------------------------------------*/

int c_register_timer(char *desc, int type)
{
   f_int fkey;
   f_int ftype;

   ftype = (f_int) type; 
   F77_NAME(register_timer, REGISTER_TIMER)(desc, &ftype, &fkey); 
   return ((int) fkey);
}

void c_update_timer(int key)
{
   f_int fkey;
   
   fkey = key;
   F77_NAME(update_timer, UPDATE_TIMER) (&fkey); 
}

void c_timer_start(int key)
{
   f_int fkey;   

   fkey = key;
   F77_NAME(timer_start, TIMER_START) (&fkey);  
}

void c_stop_timers()
{
   F77_NAME(stop_timers, STOP_TIMERS)(); 
}

void c_print_timer(char *str, int key)
{
   f_int fkey;

   fkey = key;
   F77_NAME(print_timer, PRINT_TIMER)(str, &fkey); 
}

f_int F77_NAME(get_cputime, GET_CPUTIME) ()
{
   f_int ret;
   
   ret = (f_int) clock(); 
   return ret; 
}

double F77_NAME(cpu_diff, CPU_DIFF) (f_int *ftime)
{
   /* Returns the difference (in seconds) between a previously measured 
      CPU time (argument ftime), and the current CPU time. */

   clock_t tmark, tv;
   int i;
   long max_ticks;
   clock_t diff;
   double fdiff;

   tmark = *ftime;
   tv = clock();            /* Get current clock time */
   if ( (tv < 0 && tmark < 0) ||
        (tv > 0 && tmark > 0) )
   {
      /* Signs are the same, no compensation is needed. */
      diff = tv - tmark;
   }
   else
   {
      /* Determine max. value before clock register overflows. */

      i = sizeof(clock_t);
      max_ticks = (1 << (i*8))-1;
      if (max_ticks < 0) max_ticks = (1 << (i*8-1))-1;

      /* Signs differ.  Use the following formula to compute the diff:
         diff = tv - max_ticks + max_ticks - tmark */

      diff = (tv - max_ticks) + (max_ticks - tmark);
   }

   /* Convert the difference to seconds. */
   fdiff = (double) diff / CLOCKS_PER_SEC;
   return fdiff;
}
