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
      subroutine process_server_delete_message(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_delete message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer i, j, k, n
      integer node
      integer state
      integer array, iblock, istart, nsearch

      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then
         array = server_msg(c_msg_array, node) 

c---------------------------------------------------------------------------
c   Search the server_table for entries pertaining to this array.
c---------------------------------------------------------------------------

          do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               istart = served_array_entry(i)
               nsearch = served_numblocks(i)
               go to 50
            endif
         enddo

c-----------------------------------------------------------------------------
c   There is no data for this array on the server.  Probably this is a tiny 
c   job run on more servers than necessary.  We have nothing to do.
c-----------------------------------------------------------------------------

         istart = 0
         nsearch = 0

c         print *,'Error: Cannot find array ',array,
c     *       ' in served array table'
c         call server_abort_job(server_table, nserver_table)
   50    continue

         do i = istart, istart + nsearch - 1
            if (server_table(c_server_diskloc,i) .gt. 0) then

c---------------------------------------------------------------------------
c   Negate the diskloc.  This indicates that the block is "up for grabs".
c---------------------------------------------------------------------------

               server_table(c_server_diskloc,i) = 
     *              -server_table(c_server_diskloc,i)
            endif

c---------------------------------------------------------------------------
c   Remove any references to the data in memory.
c---------------------------------------------------------------------------

            iblock = server_table(c_server_memloc,i)
            if (iblock .gt. 0) then
               if (and(server_table(c_server_flags,i),
     *                 server_dirty_flag) .ne. 0) then
                 call mark_block_clean(iblock, server_table,
     *                                  nserver_table)
               endif

               server_table(c_server_memloc,i) = 0
               server_table_ptr(iblock) = 0
            endif

c---------------------------------------------------------------------------
c   Clear the server flags for the block.
c---------------------------------------------------------------------------

            server_table(c_server_flags,i) = 0
         enddo

         server_msg(c_msg_state,node) = null_state
      endif

      return
      end
