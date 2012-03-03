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
      integer function get_symbolic_constant(ival)
c---------------------------------------------------------------------------
c   Returns the appropriate symbolic constant from the "ival" argument.
c---------------------------------------------------------------------------
      implicit none
      include 'symbolic_constants.h'
      include 'trace.h'

      integer ival
      integer value

      if (ival .ge. 0 .and. ival .lt. nsymbolic_constants) then
         value = symbolic_constant_table(ival+1)
      else
         print *,'Error: Invalid value for symbolic constant ',
     *      ival,' line ',current_line
         call abort_job()
      endif
          

      get_symbolic_constant = value
      return 
      end

      
