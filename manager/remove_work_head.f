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
      subroutine remove_work_head(node)
c---------------------------------------------------------------------------
c   Removes the head node off the server's work list, and returns its
c   node index in "node".
c
c   If the list is empty, a negative value is returned.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      integer node

      if (server_work_list_head .lt. 1) then
         node = -1
      else
         node = server_work_list_head
         server_work_list_head = server_msg(c_server_list_ptr,node) 
         server_msg(c_server_list_ptr,node) = 0

         if (server_work_list_head .eq. 0) 
     *       server_work_list_tail = 0
      endif

      return
      end

