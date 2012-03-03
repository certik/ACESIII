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
         subroutine  ofh__get_free_unit
     +
     +                    ( unitid )
     +
c------------------------------------------------------------------------
c  operation   : ofh__get_free_unit
c  module      : output files handler
c  module-id   : ofh
c  description : This operation returns the unit number of the next free
c                I/O unit. If there is no more free I/O unit available,
c                processing is aborted. This routine sets those units
c                to 'used' which have been opened in the meantime by
c                other routines in the program.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include    'output_files.h'

         logical    unit_exists

         integer    i
         integer    unitid
c
c
c------------------------------------------------------------------------
c
c
c             ...search for possible 'used' I/O units.
c
c
         do i = 2,MAX_OUTPUT_FILES
            if (file_handle (i).eq.free) then
                inquire  (unit   = i,
     +                    exist  = unit_exists)
                if (unit_exists) then
                    file_handle (i) = used
                end if
            end if
         end do
c
c
c             ...search for next free I/O unit.
c
c
         do i = 2,MAX_OUTPUT_FILES
            if (file_handle (i).eq.free) then
                unitid = i
                file_handle (i) = used_by_ofh
                return
            end if
         end do
c
c
c             ...no more I/O units available within number range.
c
c
         write (*,*) ' No more I/O units available! '
         write (*,*) ' ofh__get_free_unit '
         write (1,*) ' No more I/O units available! '
         write (1,*) ' ofh__get_free_unit '

         call  ofh__stop
c
c
c             ...ready!
c
c
         end
