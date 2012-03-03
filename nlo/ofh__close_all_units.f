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
         subroutine  ofh__close_all_units
c------------------------------------------------------------------------
c  operation   : ofh__close_all_units
c  module      : output files handler
c  module-id   : ofh
c  description : This operation terminates the connections between all
c                previously opened files and their associated units.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         include    'output_files.h'

         integer    i
c
c
c------------------------------------------------------------------------
c
c
c             ...loop over all possible I/O units and close them.
c
c
         do i = 2,MAX_OUTPUT_FILES
            call  ofh__close_unit (i)
         end do
c
c
c             ...ready!
c
c
         return
         end
