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
      subroutine check_work_stack_for_delay(flag, delay_flag,
     *                             server_msg_count)
      implicit none
      include 'server.h'

      integer i
      integer delay_flag, server_msg_count
      logical flag

      if (delay_flag .eq. 0 .or. flag .or. 
     *              barrier_in_progress) return

c--------------------------------------------------------------------------
c   We know delay_flag = 1, flag = .false. (i. e. no more probed messages
c   exist), and we are not in barrier processing mode.  
c
c   If the only messages remaining in our message nodes are requests, 
c   then we must go ahead and clear the delay_flag, otherwise the job 
c   may hang because the worker is blocked waiting for the server to 
c   clear the delay_flag and the server is blocked waiting for the request
c   to complete.
c--------------------------------------------------------------------------

      do i = 1, nmessage_buffers
         if (server_msg(c_msg_state,i) .ne. null_state .and. 
     *       server_msg(c_msg_type,i) .ne. server_request_msgtype)
     *           return
      enddo

c-------------------------------------------------------------------------
c   All active messages are requests.  Clear the flags.
c-------------------------------------------------------------------------

c      print *,'SERVER Clear delay due to requests only...'
      delay_flag = 0
      server_msg_count = 0
      return
      end


