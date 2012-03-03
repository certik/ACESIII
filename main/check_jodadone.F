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
      subroutine check_jodadone(jodadone)
c--------------------------------------------------------------------------
c   Executes the "xjoda" ACES module, and retrieves the "jodadone"
c   flag off the JOBARC file. (Only for process 0). 
c--------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'

      integer jodadone

      if (me .eq. 0) then
         call joda_post_processing(jodadone)
      endif

      return
      end
