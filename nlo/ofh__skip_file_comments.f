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
         subroutine  ofh__skip_file_comments
     +
     +                    ( unitid,comsign )
     +
c------------------------------------------------------------------------
c  operation   : ofh__skip_file_comments
c  module      : output files handler
c  module-id   : ofh
c  description : This operation skips those lines of a formatted file
c                labeled by a one character comment sign. If the file
c                as specified by its unit id is not of the type
c                sequential formatted, processing will be aborted.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include    'output_files.h'

         logical    is_open

         character*10  access_type
         character*1   char,comsign
         character*11  format_type

         integer    unitid
c
c
c------------------------------------------------------------------------
c
c
c             ...check unit number and print message/stop if necessary.
c
c
         if (unitid.lt.2  .or.  unitid.gt.MAX_OUTPUT_FILES) then
             write (*,*) ' Unit number out of range in routine: '
             write (*,*) ' ofh__skip_file_comments '
             write (1,*) ' Unit number out of range in routine: '
             write (1,*) ' ofh__skip_file_comments '
             call  ofh__stop
         end if
c
c
c             ...inquire the access type, format type and open status
c                of the unit.
c
c
         inquire  (unit   = unitid,
     +             access = access_type,
     +             form   = format_type,
     +             opened = is_open)
c
c
c             ...abort calling program, if file is not of the right
c                type.
c
c
         if ( format_type.ne.'formatted' .or.
     +        access_type.ne.'sequential' ) then
              write (*,*) ' Wrong file type! '
              write (*,*) ' access,form = ',access_type,format_type
              write (*,*) ' ofh__skip_file_comments '
              write (1,*) ' Wrong file type! '
              write (1,*) ' access,form = ',access_type,format_type
              write (1,*) ' ofh__skip_file_comments '
              call  ofh__stop
         end if
c
c
c             ...proceed if unit is open and file is being used.
c
c
         if ( is_open .and. file_handle (unitid).eq.used ) then

 1000         read (unitid,9000) char
 9000         format (a1)
              if ( char.ne.comsign ) then
                   call  ofh__backspace_unit (unitid)
                   return
              end if
              goto 1000

         else
              write (*,*) ' Unit not connected in routine: '
              write (*,*) ' ofh__skip_file_comments '
              write (1,*) ' Unit not connected in routine: '
              write (1,*) ' ofh__skip_file_comments '
              call  ofh__stop
         end if
c
c
c             ...ready!
c
c
         return
         end
