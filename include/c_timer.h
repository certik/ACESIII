#include <stdio.h>

/*-------------------------------------------------------------------------
 
   Function prototypes for the c interface of the Fortran-based timer 
   subsystem.

---------------------------------------------------------------------------*/

#define CPU_TIMER 1
#define ELAPSED_TIME_TIMER 2
#define SYSTEM_TIME_TIMER 3

int  c_register_timer(char *desc, int type);   /* Creates a timer with description */
void c_update_timer(int key);  /*  Updates a timer from the last start point */
void c_timer_start(int key);   /* Starts a timer */

