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
         subroutine  ofh__stop
c------------------------------------------------------------------------
c  operation   : ofh__stop
c  module      : output files handler
c  module-id   : ofh
c  description : This operation does all necessary operations to stop
c                carefully. All open units are closed.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include    'output_files.h'
c
c
c------------------------------------------------------------------------
c
c
c             ...finish the output files handler and stop the program.
c
c
         write (*,*) ' program stopped! '
         write (1,*) ' program stopped! '

         call  ofh__finish

         stop
c
c
c             ...ready!
c
c
         end
