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
      subroutine check_preparesum_work(node, ierr)
c-----------------------------------------------------------------------------
c   Determines if a request message may be safely started in order to insure
c   program correctness.
c
c   If there are other message nodes working on the current data block, and
c   those nodes are either prepare or preparesum messages, then this routine
c   returns a non-zero value for ierr, indicating the request may not begin. 
c   Otherwise, a 0 is returned in ierr, and the request may be safely started.
c-----------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'
      integer node, ierr
      integer i, nind, my_seqno, my_array, ptr

      my_seqno = server_msg(c_msg_seqno,node)
      my_array = server_msg(c_msg_array,node)
      nind     = server_msg(c_msg_nind,node)

c--------------------------------------------------------------------------
c   Search the list of message nodes.
c--------------------------------------------------------------------------

      ptr = server_work_list_head
  100 continue

c--------------------------------------------------------------------------
c   Is the message node a prepare or preparesum?
c---------------------------------------------------------------------------

      if (server_msg(c_msg_type,ptr) .ne. server_prepare_msgtype .and.
     *    server_msg(c_msg_type,ptr) .ne. 
     *                   server_prepare_increment) go to 200

c--------------------------------------------------------------------------
c   Is this node's sequence number earlier than the request node's seqno?
c--------------------------------------------------------------------------

      if (server_msg(c_msg_seqno,ptr) .gt. my_seqno) go to 200

c---------------------------------------------------------------------------
c   Do the array handles match?
c---------------------------------------------------------------------------

      if (server_msg(c_msg_array,ptr) .ne. my_array) go to 200

c--------------------------------------------------------------------------
c   Do all the segment ranges match?
c--------------------------------------------------------------------------

      do i = 1, nind
         if (server_msg(c_msg_bsegs+i-1,node) .ne. 
     *       server_msg(c_msg_bsegs+i-1,ptr)) go to 200
         if (server_msg(c_msg_esegs+i-1,node) .ne. 
     *       server_msg(c_msg_esegs+i-1,ptr)) go to 200
      enddo

c--------------------------------------------------------------------------
c   We have found a node which is either a prepare or preparesum on the same
c   data as that of the request node, and which has an earlier sequence.
c   Therefore, we may not start the request until this node has completed.
c   We must return ierr = 1.
c--------------------------------------------------------------------------

      ierr = 1
      print *,'Server ',me,'Unsafe to start request for node ',
     *     node,' until node ',
     *     ptr,' has completed.'
      print *,'*** DATA FOR NODE ',ptr,' *** msg_type: ',
     *   server_msg(c_msg_type,ptr),' seqno ',
     *    server_msg(c_msg_seqno,ptr),
     *   ' array ',server_msg(c_msg_array,ptr),' segs ',
     *   (server_msg(c_msg_bsegs+i-1,ptr),
     *    server_msg(c_msg_esegs+i-1,ptr),
     *   i = 1,4),' source ',server_msg(c_msg_source,ptr),' state ',
     *   server_msg(c_msg_state,ptr),' cause ',
     *   server_msg(c_msg_cause,ptr)
      return

  200 continue

c--------------------------------------------------------------------------
c   Point to the next entry in the linked list.
c--------------------------------------------------------------------------

      ptr = server_msg(c_server_list_ptr, ptr)
      if (ptr .ne. 0) go to 100

      ierr = 0
      return
      end
