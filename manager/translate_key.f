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
      subroutine translate_key(current_line, key)
c-----------------------------------------------------------------------
c   Translates a SIAL line number into a table of keys used for statistical
c   data collection.
c-------------------------------------------------------------------------
      implicit none
      include 'server_stat.h'

      integer i, current_line, key

c--------------------------------------------------------------------------
c  Searching table of line numbers for a match to the current line number.   
c--------------------------------------------------------------------------

      do i = 1, mx_stat_keys
         if (current_line .eq. lineno(i)) then
            key = i
            return
         endif

c--------------------------------------------------------------------------
c   If lineno(i) is 0, we have reached the end of the line numbers that
c   have already been entered in the table.  Create a new entry and return.
c--------------------------------------------------------------------------

         if (lineno(i) .eq. 0) then
            lineno(i) = current_line
            key = i
            return 
         endif
      enddo

c---------------------------------------------------------------------------
c   No more entries are available.
c---------------------------------------------------------------------------

      key = 0
      return
      end
