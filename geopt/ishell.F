C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.

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

