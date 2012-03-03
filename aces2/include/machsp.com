#ifndef _MACHSP_COM_
#define _MACHSP_COM_

#if defined(__fortran)

#if defined(__fortran77)
c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end

#else /* defined(__fortran9x) */
integer ::      iintln, ifltln, iintfp, ialone, ibitwd
common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
save   /machsp/

#endif /* defined(__fortran77) */

#else /* defined( __ansi_c || __cplusplus ) */

#   ifdef __cplusplus
       extern "C" {
#   endif

#include "f77_name.h"
#include "f_types.h"

#define f_machsp F77_CB_NAME(machsp,MACHSP)
struct { f_int iintln, ifltln, iintfp, ialone, ibitwd; } f_machsp;

#   ifdef __cplusplus
       }
#   endif

#endif /* defined(__fortran) */

#endif /* _MACHSP_COM_ */
