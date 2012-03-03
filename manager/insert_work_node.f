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
      subroutine insert_work_node(node)
c--------------------------------------------------------------------------
c   Inserts a node at the tail of the server's work list.
c--------------------------------------------------------------------------
      implicit none
      include 'server.h'

      integer node
      integer ptr

      if (server_work_list_tail .eq. 0) then
         server_work_list_head = node
      else
         server_msg(c_server_list_ptr,server_work_list_tail) = node
      endif

      server_work_list_tail = node
      server_msg(c_server_list_ptr,node) = 0
      return
      end
