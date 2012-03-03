#ifndef _ICORE_COM_
#define _ICORE_COM_

#if defined(__fortran)

#if defined(__fortran77)
c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end

#else /* defined(__fortran9x) */
integer :: icore(1)
common / / icore

#endif /* defined(__fortran77) */

#else /* C(++) code */
CRASH HARD!
Blank common blocks can TECHNICALLY be done, but their use is HIGHLY
ill-advised.

#endif /* defined(__fortran) */

#endif /* _ICORE_COM_ */

#ifdef _SBCORE_COM_
#error "icore.com is incompatible with sbcore.com"
#endif

