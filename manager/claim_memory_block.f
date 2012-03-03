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
      subroutine claim_memory_block(node, server_table,
     *                          nserver_table, restore_data_flag)
c---------------------------------------------------------------------------
c   Finds an available block to use for an operation.  On return, the
c   index of the block is in the nodes' c_msg_memptr field, and the 
c   next processing state is in the "c_msg_state" field.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'

      integer node, nserver_table, next_state
      integer server_table(lserver_table_entry, nserver_table)
      logical restore_data_flag

      integer iblock, ptr, ptr2, block_flag, diskloc
      integer i, busy_node

c---------------------------------------------------------------------------
c   Look for the block we need.  Is it currently in memory?
c---------------------------------------------------------------------------

      ptr = server_msg(c_msg_stptr,node)
      iblock = -1
      
c---------------------------------------------------------------------------
c   If the data is in a busy state, then some other node is working on it.
c   We must wait until it becomes free.
c---------------------------------------------------------------------------

      if (and(server_table(c_server_flags,ptr),
     *     server_busy_flag) .ne. 0) then

c--------------------------------------------------------------------------
c   If a different node has busied the block, we must wait until they are
c   through with it.
c--------------------------------------------------------------------------

         if (server_table(c_server_busy_node,ptr) .ne. node) then
            busy_node = server_table(c_server_busy_node,ptr)
c            print *,'Server ',me,'node ',node,
c     *                 ' RETRY LATER BUSY due to ',
c     *              ' node ',busy_node,' node data = ',
c     *              (server_msg(i,busy_node),i=1,lserver_msg_entry)
 
            server_msg(c_msg_state,node) = wait_for_block_state
            server_msg(c_msg_cause,node) = busy_cause
            return   ! try again later
         endif
      endif

      iblock = server_table(c_server_memloc,ptr)
      if (iblock .eq. 0) then

c---------------------------------------------------------------------------
c   We will need a new block.  Find a clean block that is not busy.
c---------------------------------------------------------------------------
 
         call find_clean_block(iblock, server_table, nserver_table)
         if (iblock .le. 0) then

c--------------------------------------------------------------------------
c   No clean blocks are available.  Force a backup and retry later.
c--------------------------------------------------------------------------

c            print *,'   NO CLEAN BLOCKS nclean = ',nclean_blocks,
c     *          ' clean_block_ptr ',clean_block_ptr
c            if (clean_block_ptr .eq. 0 .and. nclean_blocks .gt. 0) then
c               print *,'Clean blocks inconsistency'
c               call server_abort_job(server_table, nserver_table)
c            endif

            server_msg(c_msg_state,node) = wait_for_block_state
            server_msg(c_msg_cause, node) = backup_agent_cause
            return
         else
 
c---------------------------------------------------------------------------
c   We have a free clean block in memory.  Is the data actually on disk?
c---------------------------------------------------------------------------

            diskloc = server_table(c_server_diskloc,ptr)
            if (diskloc .eq. 0 .or. .not. restore_data_flag) then

c---------------------------------------------------------------------------
c   We got a clean block to use.  Set it up for the node's data.
c---------------------------------------------------------------------------

               ptr2 = server_table_ptr(iblock)
               if (ptr2 .gt. 0) then
                  server_table(c_server_flags,ptr2) =  0
                  server_table(c_server_memloc,ptr2) = 0
                  if (server_table(c_server_busy_node,ptr2) 
     *                     .ne. 0) then
                     print *,'Server ',me,' @node ',
     *                  server_table(c_server_busy_node,ptr2),
     *                  ' should have ptr2 ',ptr2,
     *                    ' busy, but it is clean'
                     call server_abort_job(server_table, nserver_table)
                  endif
               endif

               server_table(c_server_memloc,ptr) = iblock
               server_table_ptr(iblock)          = ptr
               server_table(c_server_flags,ptr)  = server_busy_flag
               server_table(c_server_busy_node,ptr) = node
               server_msg(c_msg_state,node) = null_state
               server_msg(c_msg_memptr,node) = iblock
               server_msg(c_msg_flag,node)   = preparesum_copy_flag
               return
            else

c----------------------------------------------------------------------------
c   The data is on disk.  Do we need to restore it into the memory block?
c----------------------------------------------------------------------------

               if (restore_data_flag) then

c---------------------------------------------------------------------------
c   We got a clean block to use.  Set it up for the data we are restoring.
c---------------------------------------------------------------------------

                  ptr2 = server_table_ptr(iblock)
                  if (ptr2 .gt. 0) then
                     server_table(c_server_flags,ptr2) =  0
                     server_table(c_server_memloc,ptr2) = 0
                     if (server_table(c_server_busy_node,ptr2) 
     *                     .ne. 0) then
                        print *,'Server ',me,' #node ',
     *                  server_table(c_server_busy_node,ptr2),
     *                  ' should have ptr2 ',ptr2,
     *                    ' busy, but it is clean'
                        call server_abort_job(server_table, 
     *                       nserver_table)
                     endif
                  endif

                  server_table(c_server_memloc,ptr) = iblock
                  server_table_ptr(iblock)          = ptr
                  server_table(c_server_flags,ptr) = 
     *                                       server_busy_flag
                  server_table(c_server_busy_node,ptr) = 
     *                                       node
                  server_msg(c_msg_state,node) = wait_for_block_state
                  server_msg(c_msg_cause,node) = restore_cause
                  server_msg(c_msg_memptr,node) = iblock
c                  print *,'   NODE ',node,' WAIT FOR RESTORE iblock ',
c     *                   iblock
                  return
               endif   ! restore_data_flag
            endif      ! diskloc .eq. 0...
         endif         ! find_clean_block
      else             ! iblock .eq. 0

c---------------------------------------------------------------------------
c   The block is in memory and  available to be used.
c   Mark the block busy, save it in the message node.
c----------------------------------------------------------------------------

         server_msg(c_msg_memptr,node) = iblock
         server_table(c_server_flags,ptr) = or(server_busy_flag,
     *                          server_table(c_server_flags,ptr))
         server_table(c_server_busy_node,ptr) = node
         server_msg(c_msg_state,node) = null_state
         server_msg(c_msg_cause,node) = null_cause
c         print *,'SUCCESSFUL CLAIM FOR NODE ',node,' type ',
c     *     server_msg(c_msg_type,node),' block ',iblock
      endif
 
      return
      end
