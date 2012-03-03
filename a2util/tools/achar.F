
c This routine supplies the character*1 function 'achar' if the Fortran
c compiler does not contain it as an intrinsic. This was first seen on
c an OSF/DEC compiler.

c Recall: achar returns the ASCII character of a decimal number while
c         char returns the internal character representation of a decimal.

#ifdef _CHAR_IS_ACHAR
      character*1 function achar(num)
      implicit none
      integer num
      character*1 char
      external char
      achar = char(num)
      return
      end
#else
#include "return.f"
#endif /* _CHAR_IS_ACHAR */

