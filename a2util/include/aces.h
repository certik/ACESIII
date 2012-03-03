
#ifndef _ACES_H_
#define _ACES_H_ /* ACES DATA TYPE HANDLES AND SYSTEM-WIDE DEFINES */

#undef M_REAL
#undef M_SINGLE
#undef M_DOUBLE

#undef M_IMPLICIT
#undef M_TRACEBACK

#define F_INTEGER   0
#define F_REAL      1
#define F_COMPLEX   2
#define F_LOGICAL   3
#define F_CHARACTER 4

#ifdef __fortran

#   ifdef __fortran77

c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler

#define M_REAL double precision
#define M_DOUBLE
#define M_IMPLICIT implicit none

#include <blas.h>
#include <matx.h>
cYAU - ACES3 stuff . . . we hope - #include <aces.par>

#   else /* __fortran77 */

! INSERT FORTRAN9x CODE HERE

#   endif /* __fortran77 */

#else /* __fortran */

#define M_REAL double

#endif /* __fortran */

#endif /* _ACES_H_ */

