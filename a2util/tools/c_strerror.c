
/*
 * This routine maps the global error number into an error message and
 * prints the string to stdout and stderr.
 */

#include <errno.h>
#include <stdio.h>
#include <string.h>

#include "f77_name.h"

void
F77_NAME(c_strerror,C_STRERROR)
()
{
    char *p = strerror(errno);
    fprintf(stdout," System error: %d, %s\n",errno,p);
    fprintf(stderr," System error: %d, %s\n",errno,p);
    fflush(stdout);
    fflush(stderr);
}

