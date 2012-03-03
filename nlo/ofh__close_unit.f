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
         subroutine  ofh__close_unit 
     +
     +                    ( unitid )
     +
c------------------------------------------------------------------------
c  operation   : ofh__close_unit
c  module      : output files handler
c  module-id   : ofh
c  description : This operation terminates the connection between the
c                specified unit and the associated file. If the unit
c                number is within range mark it as free.
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

         integer    unitid
c
c
c------------------------------------------------------------------------
c
c
c             ...inquire the open status of the unit.
c
c
         inquire  (unit   = unitid,
     +             opened = is_open)
c
c
c             ...close the unit if it is open.
c
c
         if (is_open .and.
     *       file_handle(unitid) .eq. used_by_ofh) then
             close (unit = unitid)
         endif
c
c
c             ...mark the unit as free if within range.
c
c
         if (unitid.ge.2  .or.  unitid.le.MAX_OUTPUT_FILES) then
             if (file_handle(unitid) .eq. used_by_ofh)
     *           file_handle (unitid)  =  free
         end if
c
c
c             ...ready!
c
c
         return
         end
