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
         subroutine  ofh__open_file
     +
     +                    ( file_name,
     +                      file_status,
     +                      file_access,
     +                      file_format,
     +
     +                              unitid )
     +
c------------------------------------------------------------------------
c  operation   : ofh__open_file
c  module      : output files handler
c  module-id   : ofh
c  description : This operation opens a file and returns the unit number.
c                Processing is aborted, if one of the following cases
c                apply:
c
c                     1) file is already open
c                     2) file should have status old, but is
c                        not existing
c                     3) file should have status new, but is
c                        existing
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include    'output_files.h'

         character*(*)      file_name
         character*(*)      file_status
         character*(*)      file_access
         character*(*)      file_format

         logical            file_exists
         logical            file_open

         integer            record_length
         integer            unitid
c
c
c------------------------------------------------------------------------
c
c
c             ...check for proper file attributes.
c
c
         if (file_status .ne. 'old'     .or.
     +       file_status .ne. 'new'     .or.
     +       file_status .ne. 'scratch' .or.
     +       file_status .ne. 'unknown'      ) then
             write (*,*) ' Unknown file status! '
             write (*,*) ' File status = ',file_status
             write (*,*) ' Cannot open file! '
             write (*,*) ' ofh__open_file '
             write (1,*) ' Unknown file status! '
             write (1,*) ' File status = ',file_status
             write (1,*) ' Cannot open file! '
             write (1,*) ' ofh__open_file '
             call  ofh__stop
         end if

         if (file_access .ne. 'direct'     .or.
     +       file_access .ne. 'sequential'      ) then
             write (*,*) ' Unknown file access method! '
             write (*,*) ' File access = ',file_access
             write (*,*) ' Cannot open file! '
             write (*,*) ' ofh__open_file '
             write (1,*) ' Unknown file access method! '
             write (1,*) ' File access = ',file_access
             write (1,*) ' Cannot open file! '
             write (1,*) ' ofh__open_file '
             call  ofh__stop
         end if

         if (file_format .ne. 'formatted'   .or.
     +       file_format .ne. 'unformatted'      ) then
             write (*,*) ' Unknown file format! '
             write (*,*) ' File format = ',file_format
             write (*,*) ' Cannot open file! '
             write (*,*) ' ofh__open_file '
             write (1,*) ' Unknown file format! '
             write (1,*) ' File format = ',file_format
             write (1,*) ' Cannot open file! '
             write (1,*) ' ofh__open_file '
             call  ofh__stop
         end if
c
c
c             ...inquire the open and existence status of the file
c                associated with the file name.
c
c
         inquire  (file   = file_name,
     +             opened = file_open,
     +             exist  = file_exists)
c
c
c             ...file already open.
c
c
         if (file_open) then
             write (*,*) ' The file: ',file_name,' is already open! '
             write (*,*) ' ofh__open_file '
             write (1,*) ' The file: ',file_name,' is already open! '
             write (1,*) ' ofh__open_file '
             call  ofh__stop
         end if
c
c
c             ...check file existence if status = old or new.
c
c
         if (.not.file_exists  .and.  file_status.eq.'old') then
             write (*,*) ' The file: ',file_name,' does not exist! '
             write (*,*) ' ofh__open_file '
             write (1,*) ' The file: ',file_name,' does not exist! '
             write (1,*) ' ofh__open_file '
             call  ofh__stop
         end if

         if (file_exists  .and.  file_status.eq.'new') then
             write (*,*) ' The file: ',file_name,' already exists! '
             write (*,*) ' ofh__open_file '
             write (1,*) ' The file: ',file_name,' already exists! '
             write (1,*) ' ofh__open_file '
             call  ofh__stop
         end if
c
c
c             ...if scratch file is wanted set the record length.
c
c
         if (file_access .eq. 'direct') then
             record_length = 4096
         end if
c
c
c             ...get next free unit number.
c
c
         call  ofh__get_free_unit (unitid)
c
c
c             ...open I/O unit.
c
c
         if (file_status .eq. 'scratch') then

             open  (unit   = unitid,
     +              access = file_access,
     +              form   = file_format,
     +              status = 'scratch')

         else

             if (file_access .eq. 'direct') then

                 open  (unit   = unitid,
     +                  access = 'direct',
     +                  file   = file_name,
     +                  form   = file_format,
     +                  recl   = record_length,
     +                  status = file_status)

             else

                 open  (unit   = unitid,
     +                  access = file_access,
     +                  file   = file_name,
     +                  form   = file_format,
     +                  status = file_status)

             end if

         end if
c
c
c             ...ready!
c
c
         return
         end
