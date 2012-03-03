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
         subroutine  ofh__open_logfile
c------------------------------------------------------------------------
c  operation   : ofh__open_logfile
c  module      : output files handler
c  module-id   : ofh
c  description : This operation opens the logfile for sequential access
c                and formatted output and connects it to unit 1. If the
c                logfile is already open, the routine aborts the
c                program.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include    'output_files.h'

         logical    file_open
c
c
c------------------------------------------------------------------------
c
c
c             ...inquire the open status of the logfile.
c
c
         inquire  (unit   = 1,
     +             opened = file_open )
c
c
c             ...logfile already open.
c
c
         if (file_open) then
             write (*,*) ' Cannot connect Logfile to unit 1 ! '
             write (*,*) ' Program will abort! '
             write (*,*) ' ofh__open_logfile '
             call  ofh__stop
         end if
c
c
c             ...open logfile unit and write header.
c
c
         open  (unit   = 1,
     +          access = 'sequential',
     +          file   = logfile_name,
     +          form   = 'formatted',
     +          status = 'unknown')

         write (1,9000) ' ACES3 LOGFILE '
 9000    format (///,a15,///)

         logfile_exists = .true.
c
c
c             ...ready!
c
c
         return
         end
