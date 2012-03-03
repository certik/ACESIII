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
         subroutine  ofh__backspace_unit
     +
     +                    ( unitid )
     +
c------------------------------------------------------------------------
c  operation   : ofh__backspace_unit
c  module      : output files handler
c  module-id   : ofh
c  description : This operation performs a backspace command on the
c                unit specified by its unit id.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include       'output_files.h'

         logical       is_open

         character*10  access_type

         integer       unitid
c
c
c------------------------------------------------------------------------
c
c
c             ...inquire the access type and open status of the unit.
c
c
         inquire (unit   = unitid,
     +            access = access_type,
     +            opened = is_open)
c
c
c             ...abort calling program, if file is not of the right
c                type.
c
c
         if (access_type.ne.'sequential') then
             write (*,*) ' Wrong file type! '
             write (*,*) ' access type = ',access_type
             write (*,*) ' ofh__backspace_unit '
             write (1,*) ' Wrong file type! '
             write (1,*) ' access type = ',access_type
             write (1,*) ' ofh__backspace_unit '
             call  ofh__stop
         end if
c
c
c             ...proceed if unit is open and file is being used.
c
c
         if (is_open .and. file_handle (unitid).eq.used) then
             backspace (unit = unitid)
         else
             write (*,*) ' Unit not connected! '
             write (*,*) ' ofh__backspace_unit '
             write (1,*) ' Unit not connected! '
             write (1,*) ' ofh__backspace_unit '
             call  ofh__stop
         end if
c
c
c             ...ready!
c
c
         return
         end
