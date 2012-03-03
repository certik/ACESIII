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
         subroutine  ofh__finish
c------------------------------------------------------------------------
c  operation   : ofh__finish
c  module      : output files handler
c  module-id   : ofh
c  description : This operation finishes the output files handler.
c                The following is done in this routine:
c
c                       1) the logfile is closed.
c                       2) all current open units are closed.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
         include 'output_files.h'
c
c
c------------------------------------------------------------------------
c
c
c             ...close the logfile and all units.
c
c
         if (logfile_exists) then
             call  ofh__close_logfile
         end if

         call  ofh__close_all_units
c
c
c             ...print information on screen.
c
c
         write (*,*) ' Logfile and all units closed! '
c
c
c             ...ready!
c
c
         return
         end
