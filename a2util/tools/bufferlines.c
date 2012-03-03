
/*
 * This routine changes the size of the output buffer.
 */

#include <stdio.h>

/*
   YAU - Someone check this. `man setvbuf` says the NULL buffer is non-portable.
*/

/*
void bufferlines ()
{ setvbuf(stdout,NULL,_IOLBF,1024); return; }

void bufferlines_()
{ setvbuf(stdout,NULL,_IOLBF,1024); return; }

void BUFFERLINES ()
{ setvbuf(stdout,NULL,_IOLBF,1024); return; }
*/

void bufferlines ()
{ setlinebuf(stdout); return; }

void bufferlines_()
{ setlinebuf(stdout); return; }

void BUFFERLINES ()
{ setlinebuf(stdout); return; }

