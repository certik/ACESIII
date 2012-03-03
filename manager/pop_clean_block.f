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
      subroutine pop_clean_block(iblock, server_table, nserver_table)
c------------------------------------------------------------------------
c   Pops a block off the clean block stack.  Returns iblock = -1 if
c   the stack is empty.
c-------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'

      integer iblock, nserver_table
      integer server_table(lserver_table_entry,nserver_table)

      integer i, j, jblock, flagval, ptr

c--------------------------------------------------------------------------
c   Search backwards through the clean block "stack" until we find one that is
c   not busy.  
c--------------------------------------------------------------------------

      iblock = -1

      if (clean_block_ptr .lt. 1) go to 200 

      flagval = or(server_busy_flag, server_dirty_flag)

      do i = clean_block_ptr, 1, -1
         jblock = clean_blocks(i)
         ptr    = server_table_ptr(jblock)
         if (ptr .eq. 0 .or.
     *      (ptr .gt. 0 .and.
     *       and(flagval, server_table(c_server_flags,ptr)) .eq.
     *                           0)) then
            iblock = jblock
            j      = i
            go to 100
         endif
      enddo

c---------------------------------------------------------------------------
c   There are no clean blocks that are not busy.  
c---------------------------------------------------------------------------

      go to 200     

  100 continue

c----------------------------------------------------------------------------
c   Remove iblock from the list of clean blocks.
c----------------------------------------------------------------------------

      do i = j+1, clean_block_ptr
         clean_blocks(i-1) = clean_blocks(i) 
      enddo

      clean_blocks(clean_block_ptr) = 0
      clean_block_ptr = clean_block_ptr - 1
      nclean_blocks = nclean_blocks - 1

  200 continue
      return
      end
