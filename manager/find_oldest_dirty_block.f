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
      subroutine find_oldest_dirty_block(iblock, server_table,
     *                  nserver_table)
c---------------------------------------------------------------------------
c   Finds the oldest dirty block in memory that is also not busy.
c   If none exists, returns iblock = -1.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      integer iblock
      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
 
      integer ptr

      if (dirty_list_head .eq. 0) then
         iblock = -1
         return
      endif

      iblock = dirty_list_head
  100 continue
      ptr = server_table_ptr(iblock)
      if (and(server_table(c_server_flags, ptr), server_busy_flag) .eq.
     *                 0) return
 
      iblock = dirty_list_ptr(iblock)
      if (iblock .ne. 0) go to 100

      iblock = -1
      return
      end
