
c This routine supplies the integer function 'iachar' if the Fortran
c compiler does not contain it as an intrinsic.

c Recall: iachar returns the decimal number of an ASCII character while
c         ichar returns the internal decimal of a character representation.

#ifdef _ICHAR_IS_IACHAR
      integer function iachar(cz)
      implicit none
      character*1 cz
      integer ichar
      external ichar
      iachar = ichar(cz)
      return
      end
#else
#include "return.f"
#endif /* _ICHAR_IS_IACHAR */

