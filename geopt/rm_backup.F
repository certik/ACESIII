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

c This routine removes the CURRENT backup directory.

#include "flags.h"

      subroutine rm_backup
      implicit none

      character*80 szSaveDir
      integer i, iTmp
#include "jodaflags.com"

      if (iflags(h_IFLAGS_restart).ne.0) then
         call gfname('SAVEDIR',szSaveDir,i)
         print *, 'Removing ',szSaveDir(1:i),'/CURRENT . . .'
c         iTmp=ishell('/bin/rm -rf '//szSaveDir(1:i)//'/CURRENT')
         if (iTmp.eq.0) then
            print *, 'The backup directory was successfully removed.'
         else
            print *, '@RM_BACKUP: The directory ',szSaveDir(1:i),
     &               '/CURRENT could not be removed.'
            print *, '          Future calculations may be ',
     &               'contaminated.'
         end if
c     end if (restart)
      end if

      return
c     end subroutine rm_backup
      end

