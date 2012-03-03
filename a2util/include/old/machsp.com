#ifndef _MACHSP_COM_
#define _MACHSP_COM_

#if defined(__fortran)

#if defined(__fortran77)
c machsp.com : begin
c iintln := byte-length of an integer
c ifltln := byte-length of a double precision float
c iintfp := number of integers in a double precision float
c ialone := 2**(8*iintln/4)-1
c ibitwd := 8*iintln/4
      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
c machsp.com : end

#else /* defined(__fortran9x) */
integer ::      iintln, ifltln, iintfp, ialone, ibitwd
common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd

#endif /* defined(__fortran77) */

#else /* defined( __ansi_c || __cplusplus ) */

#   ifdef __cplusplus
       extern "C" {
#   endif

#include "f_types.h"

struct struct_machsp
{
    f_int iintln, ifltln, iintfp, ialone, ibitwd;
};
#ifdef C_SUFFIX
   struct struct_machsp machsp_;
#  define MACHSP machsp_
#else
   struct struct_machsp machsp;
#  define MACHSP machsp
#endif /* C_SUFFIX */

#   ifdef __cplusplus
       }
#   endif

#endif /* defined(__fortran) */

#endif /* _MACHSP_COM_ */
