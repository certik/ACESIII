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
      subroutine make_free_block(node, server_table, nserver_table)
c-----------------------------------------------------------------------------
c   Produces a free block, doing a backup of a dirty block to disk if 
c   necessary.  If successful, the node's "c_msg_state" field is set to
c   "null_state", and the "c_msg_memptr" field points to the block.
c-----------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_stat.h'
      include 'dbugcom.h'

      integer nserver_table, node
      integer server_table(lserver_table_entry,nserver_table)
      integer iblock

c---------------------------------------------------------------------------
c   First, see if a clean block is already available.
c---------------------------------------------------------------------------

      call find_clean_block(iblock, server_table, nserver_table)
      if (iblock .le. 0) then

c----------------------------------------------------------------------------
c   All clean blocks are busy.  We must find a dirty block and back it up 
c   to disk.
c---------------------------------------------------------------------------

         call do_backup(iblock, server_table, nserver_table)
      endif

      if (iblock .gt. 0) then
         server_msg(c_msg_state,node)  = null_state
         server_msg(c_msg_memptr,node) = iblock
      endif
      return
      end
