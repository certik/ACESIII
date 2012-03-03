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
         subroutine  ofh__start
c------------------------------------------------------------------------
c  operation   : ofh__start
c  module      : output files handler
c  module-id   : ofh
c  description : This operation starts the output files handler.
c                The following is done in this routine:
c
c                       1) the logfile is opened.
c                       2) all unit numbers (identifiers) currently
c                          not in use are set free.
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
c
c
c------------------------------------------------------------------------
c
c
c             ...set initial values.
c
c
         logfile_exists = .false.

         do i = 2,MAX_OUTPUT_FILES
            file_handle (i)  =  free
         end do
c
c
c             ...open logfile.
c
c
         call  ofh__open_logfile
c
c
c             ...set all i/o units to free with exception of:
c
c                     unit 1  =  logfile
c                     unit 5  =  standard input unit (system)
c                     unit 6  =  standard output unit (system)
c                     unit x  =  all units which are currently in use
c
c
         file_handle (1) = used
         file_handle (5) = used
         file_handle (6) = used

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
c             ...print information on screen.
c
c
         write (*,*) ' Logfile opened and all I/O units ready! '
c
c
c             ...ready!
c
c
         return
         end
