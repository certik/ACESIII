
#if defined(_UNICOS) || defined(_IRIX)
c Use system intrinsic
c (use -craylibs flag on SGI link)
#include "return.f"
#else
      integer function ishell(string)
      character*(*) string
      character*1 achar
      integer c_system
      external c_system
      character*80 cmd
      if ( len(string) .lt. 80 ) then
         cmd = string//achar(0)
         ishell = c_system(cmd)
      else
         print *, '@ISHELL: The command buffer is only 80 characters.'
         print *, '         exiting while trying to execute: ',
     &            '(string between arrows)'
         print '(3a)', '-->',string,'<--'
         call c_exit(1)
      end if
      return
      end

#endif /* defined(_UNICOS) || defined(_IRIX) */

