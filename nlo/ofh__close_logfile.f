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
         subroutine  ofh__close_logfile
c------------------------------------------------------------------------
c  operation   : ofh__close_logfile
c  module      : output files handler
c  module-id   : ofh
c  description : This operation terminates the connections between the
c                unit 1 and the logfile. the contents of the logfile
c                will be saved.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         logical    is_open
c
c
c------------------------------------------------------------------------
c
c
c             ...inquire the open status of the logfile unit.
c
c
         inquire  (unit   = 1,
     +             opened = is_open)
c
c
c             ...close the logfile unit if it is open.
c
c
         if (is_open) then
             close (unit = 1)
         endif
c
c
c             ...ready!
c
c
         return
         end
