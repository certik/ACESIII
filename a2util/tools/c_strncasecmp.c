
/*
 * This routine compares two strings up to len bytes ignoring case.
 */

#include <strings.h>

#include "f_types.h"

f_int c_strncasecmp (const char * s1, const char * s2, f_int * len)
{ return (f_int)strncasecmp(s1,s2,(size_t)*len); }

f_int c_strncasecmp_(const char * s1, const char * s2, f_int * len)
{ return (f_int)strncasecmp(s1,s2,(size_t)*len); }

f_int C_STRNCASECMP (const char * s1, const char * s2, f_int * len)
{ return (f_int)strncasecmp(s1,s2,(size_t)*len); }

