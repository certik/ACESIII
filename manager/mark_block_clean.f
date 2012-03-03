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
      subroutine mark_block_clean(iblock, server_table,nserver_table)
c----------------------------------------------------------------------------
c   Marks a block as "clean".
c----------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'
      integer iblock, nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer ptr, flagval

      ptr = server_table_ptr(iblock)
      flagval = and(server_table(c_server_flags,ptr),
     *              server_dirty_flag)
      server_table(c_server_flags,ptr) = xor(flagval,
     *       server_table(c_server_flags,ptr))
      call push_clean_block(iblock, server_table, nserver_table)

c-----------------------------------------------------------------------------
c   Remove block from dirty list.
c-----------------------------------------------------------------------------

      if (dirty_list_head .eq. 0)  return
      if (dirty_list_head .eq. iblock) then
         dirty_list_head = dirty_list_ptr(iblock)
         if (iblock .eq. dirty_list_tail) dirty_list_tail = 0
         ndirty = ndirty - 1
         return
      endif 

      ptr = dirty_list_head
  100 continue
      if (ptr .eq. 0) then
         print *,'Error: Cannot find block ',iblock,
     *      ' in list of dirty blocks'
         call server_abort_job(server_table, nserver_table)
      endif

      if (dirty_list_ptr(ptr) .eq. iblock) then
         dirty_list_ptr(ptr) = dirty_list_ptr(iblock)
         dirty_list_ptr(iblock) = 0 
         if (iblock .eq. dirty_list_tail) dirty_list_tail = ptr
         ndirty = ndirty - 1
         return
      endif

      ptr = dirty_list_ptr(ptr)

      go to 100
      end

