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
         integer function  ofh__get_file_unit
     +
     +                          ( name )
     +
c------------------------------------------------------------------------
c  operation   : ofh__get_file_unit
c  module      : output files handler
c  module-id   : ofh
c  description : This operation gets the unit number for the file with
c                filename name.
c
c  author      : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...include files and declare variables.
c
c
         character*(*)  name
c
c
c------------------------------------------------------------------------
c
c
c             ...inquire the unit number.
c
c
         inquire  ( file   = name ,
     +              number = ofh__get_file_unit )
c
c
c             ...ready!
c
c
         return
         end
