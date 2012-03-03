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
      subroutine push_work_node(node)
c-------------------------------------------------------------------------
c   This subroutine pushes a node onto the stack of available message nodes.
c   The index of the node is specified in the "node" argument.
c
c   After the node has been pushed back onto the stack, it is available
c   for use again.
c-------------------------------------------------------------------------
      implicit none
      include 'server.h'

      integer node
     
      if (server_node_ptr .ge. mx_server_msg) then
         print *,'Error: Server message stack overflow'
         call abort_job() 
      else
         if (node .lt. 1 .or. node .gt. mx_server_msg) then
            print *,'Error: Invalid data on server stack'
            print *,'   node = ',node,' should be >= 1 and <= ',
     *               mx_server_msg
            call abort_job()
         else
            server_node_ptr = server_node_ptr + 1
            server_msg_node(server_node_ptr) = node
         endif
      endif

      return
      end

