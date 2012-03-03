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
      integer function find_free_diskloc(ifile, server_table, 
     *                                   nserver_table)
c--------------------------------------------------------------------------
c   Finds a free disk location for the correct file.  If no free blocks are
c   available already, this function returns -1.
c--------------------------------------------------------------------------
      implicit none
      include 'server.h'

      integer ifile, nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer i

      find_free_diskloc = -1

      do i = 1, nserver_table
         if (server_table(c_server_file,i) .eq. ifile) then
            if (server_table(c_server_diskloc,i) .lt. 0) then
               find_free_diskloc = -server_table(c_server_diskloc,i)
               server_table(c_server_diskloc,i) = 0
               return
            endif
         endif
      enddo

      return
      end
