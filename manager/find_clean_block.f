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
      subroutine find_clean_block(iblock, server_table, nserver_table)
c---------------------------------------------------------------------------
c   Find a block that is in a "clean" state.
c   If no clean blocks are found that are not also busy, return -1.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'

      integer iblock, nserver_table
      integer server_table(lserver_table_entry,nserver_table)

      call pop_clean_block(iblock, server_table, nserver_table)

      return
      end
