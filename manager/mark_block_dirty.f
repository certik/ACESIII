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
      subroutine mark_block_dirty(iblock, server_table,nserver_table)
c----------------------------------------------------------------------------
c   Marks a block as "dirty".
c----------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'
      integer iblock, nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer ptr
      integer jblock
      integer save_flags

      ptr = server_table_ptr(iblock)

c----------------------------------------------------------------------------
c   If the block is already dirty, remove it from the dirty list, then
c   add it to the tail of the dirty list.  Since the algorithm to make free
c   blocks works by taking the oldest dirty block, this will tend to 
c   retain blocks in the server's memory longer.
c----------------------------------------------------------------------------

      if (and(server_table(c_server_flags,ptr),
     *        server_dirty_flag) .ne. 0) then 
         call remove_from_dirty_list(iblock, 
     *                               server_table,nserver_table) 
      endif

c---------------------------------------------------------------------------
c   Now turn on the dirty flag.
c---------------------------------------------------------------------------

      server_table(c_server_flags,ptr) = or(server_dirty_flag,
     *       server_table(c_server_flags,ptr))

c-----------------------------------------------------------------------------
c   Add to the dirty block list at the tail.
c-----------------------------------------------------------------------------

      if (dirty_list_tail .eq. 0) then
         dirty_list_tail = iblock
         dirty_list_head = iblock
      else
         ptr = dirty_list_tail
         dirty_list_tail = iblock
         dirty_list_ptr(ptr) = iblock
      endif

      dirty_list_ptr(iblock) = 0
      ndirty = ndirty + 1
      return
      end
