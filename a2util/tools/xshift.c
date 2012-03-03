
/*
 * These routines right/left shift an integer i by a variable number of bits n.
 */

#include "f_types.h"

f_int rshift (f_int * i, f_int * n)
{ return (f_int)(*i >> *n); }

f_int rshift_(f_int * i, f_int * n)
{ return (f_int)(*i >> *n); }

f_int RSHIFT (f_int * i, f_int * n)
{ return (f_int)(*i >> *n); }

f_int lshift (f_int * i, f_int * n)
{ return (f_int)(*i << *n); }

f_int lshift_(f_int * i, f_int * n)
{ return (f_int)(*i << *n); }

f_int LSHIFT (f_int * i, f_int * n)
{ return (f_int)(*i << *n); }

