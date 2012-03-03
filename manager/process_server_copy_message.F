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
      subroutine process_server_copy_message(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_copy message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer i, j, ii
      integer node
      integer state
      integer match, nind
      integer source_array, target_array
      integer sentry, tentry
      integer src_bseg(mx_array_index)
      integer src_eseg(mx_array_index)
      integer tgt_bseg(mx_array_index)
      integer tgt_eseg(mx_array_index)

      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then
         source_array = server_msg(c_msg_nind, node)
         target_array = server_msg(c_msg_array, node) 

         do i= 1, nserver_table

c------------------------------------------------------------------------
c   Find the next block of the source_array in the server_table.
c------------------------------------------------------------------------

            if (server_table(c_server_array,i) .eq. 
     *                                   source_array) then

c------------------------------------------------------------------------
c   Save the block's indices.
c------------------------------------------------------------------------

               nind = server_table(c_server_nind,i)
               do j = 1, nind
                  src_bseg(j) = server_table(c_server_bsegs+j-1,i)
                  src_eseg(j) = server_table(c_server_esegs+j-1,i)
               enddo

c------------------------------------------------------------------------
c   Find the corresponding block of the target array.
c------------------------------------------------------------------------

               do ii = 1, nserver_table
                  if (server_table(c_server_array,ii) .eq.
     *                target_array) then
                    
c-------------------------------------------------------------------------
c   Check for a match to the indices of the source block.
c-------------------------------------------------------------------------

                    match = 0
                    do j = 1, nind
                       if (src_bseg(j) .eq. 
     *                    server_table(c_server_bsegs+j-1,ii) .and.
     *                     src_eseg(j) .eq. 
     *                    server_table(c_server_esegs+j-1,ii) )
     *                      match = match + 1   
                    enddo

                    if (match .eq. nind) then

c------------------------------------------------------------------------
c   We have a matching block.  We must swap the array information.
c------------------------------------------------------------------------

                       server_table(c_server_array,i) = -target_array
                       server_table(c_server_array,ii) = -source_array
                       go to 100
                    endif
                  endif 
               enddo

c-------------------------------------------------------------------------
c   We never found a matching block for the source block.  Thus it lives 
c   on another server.  We simply relabel this block as the target.
c-------------------------------------------------------------------------

               server_table(c_server_array,i) = -target_array
            endif

  100    continue
         enddo

c------------------------------------------------------------------------
c   Now change back the flagged array values to the correct array indices.
c------------------------------------------------------------------------

         do i = 1, nserver_table
c-------------------------------------------------------------------------
c   Change any blocks still labeled as target_array to source_array.
c   They are blocks that were never changed because their "mate" lives
c   on another server.
c-------------------------------------------------------------------------

            if (server_table(c_server_array,i) .eq. target_array)
     *          server_table(c_server_array,i) = source_array

c--------------------------------------------------------------------------
c   Reverse the sign on the changed blocks.
c--------------------------------------------------------------------------

            if (server_table(c_server_array,i) .eq. -target_array) 
     *          server_table(c_server_array,i) = target_array
            if (server_table(c_server_array,i) .eq. -source_array)
     *          server_table(c_server_array,i) = source_array

         enddo

c-------------------------------------------------------------------------
c   Swap the entries in the served_array_table.
c-------------------------------------------------------------------------

         sentry = 0
         tentry = 0
         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. source_array) sentry = i
            if (served_array_table(i) .eq. target_array) tentry = i
         enddo

         if (sentry .ne. 0) served_array_table(sentry) = 
     *                                 target_array
         if (tentry .ne. 0) served_array_table(tentry) = 
     *                                 source_array
  
         server_msg(c_msg_state,node) = null_state
      endif

      return
      end
