
#ifndef _JODAFLAGS_COM_
#define _JODAFLAGS_COM_

#include "flags.h"

#ifdef __fortran

#   ifdef __fortran77

c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

#ifdef MONSTER_FLAGS
      integer        ioppar(600)
      common /flags/ ioppar
      save   /flags/
#else
      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/
#endif

#   else /* __fortran77 */

  CRASH HARD!!! (since joda is not written in Fortran 90)

#   endif /* __fortran77 */

#else /* __fortran */

#   ifdef __cplusplus
       extern "C" {
#   endif

#include "f77_name.h"
#include "f_types.h"

#define f_flags F77_CB_NAME(flags,FLAGS)
#ifdef MONSTER_FLAGS
struct { f_int ioppar[600]; } f_flags;
#else
struct { f_int iflags[100], iflags2[500]; } f_flags;
#endif

#   ifdef __cplusplus
       }
#   endif

#endif /* __fortran */

#endif /* _JODAFLAGS_COM_ */

