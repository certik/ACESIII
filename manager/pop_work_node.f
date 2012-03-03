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
      subroutine pop_work_node(node)
c-------------------------------------------------------------------------
c   This subroutine pops a node off the stack of available message nodes.
c   The index of the node is returned in the "node" argument.
c
c   If no nodes are available, a negative value is returned.
c-------------------------------------------------------------------------
      implicit none
      include 'server.h'

      integer node
     
      if (server_node_ptr .lt. 1) then
         node = -1
         return
      else
         node = server_msg_node(server_node_ptr)
         server_node_ptr = server_node_ptr - 1
      endif

      return
      end

