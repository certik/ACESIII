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
         subroutine  ofh__rewind_unit
     +
     +                    ( unitid )
     +
c------------------------------------------------------------------------
c  operation   : ofh__rewind_unit
c  module      : output files handler
c  module-id   : ofh
c  description : This operation rewinds the specified unit. If the unit
c                is not connected to the unit number, processing will
c                be aborted.
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
c             ...check unit number and print message/stop if necessary.
c
c
         if ( unitid.lt.2  .or.  unitid.gt.MAX_OUTPUT_FILES ) then
              write (*,*) ' Unit number out of range in routine: '
              write (*,*) ' ofh__rewind_unit '
              write (1,*) ' Unit number out of range in routine: '
              write (1,*) ' ofh__rewind_unit '
              call  ofh__stop
         end if
c
c
c             ...inquire the open status of the unit.
c
c
         inquire  (unit   = unitid,
     +             opened = is_open)
c
c
c             ...rewind unit.
c
c
         if ( is_open  .and.  file_handle (unitid) .eq. used ) then
              rewind ( unit = unitid )
         else
              write (*,*) ' Unit not connected in routine: '
              write (*,*) ' ofh__rewind_unit '
              write (1,*) ' Unit not connected in routine: '
              write (1,*) ' ofh__rewind_unit '
              call  ofh__stop
         end if
c
c
c             ...ready!
c
c
         return
         end
