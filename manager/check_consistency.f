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
      subroutine check_consistency(call_marker,server_table, 
     *     nserver_table)
      implicit none
      include 'server.h'
      include 'parallel_info.h'

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer call_marker
      integer i, j,k,memloc
      integer iblock, ptr, node

      do i = 1, nserver_table
         ptr = i
         if (and(server_table(c_server_flags,ptr),
     *               server_busy_flag) .ne. 0) then
            node = server_table(c_server_busy_node,i)
            if (node .eq. 0) then
               print *,'Server ',me,' Error: ',call_marker,' ptr ',
     *             ptr,' is busy  but busy_node = 0'
               call server_abort_job(server_table, nserver_table)
            endif

            iblock = server_msg(c_msg_memptr,node)
            if (iblock .gt. 0) then
               if (server_table_ptr(iblock) .ne. ptr) then
                  print *,'Server ',me,' Error: ',call_marker,
     *                ' server_table_ptr inconsistent'
                  print *,'iblock, server_table_ptr(iblock) ',
     *                iblock, server_table_ptr(iblock),' should be ',
     *                ptr
                  call server_abort_job(server_table, nserver_table)
               endif

               if (server_table(c_server_memloc,ptr) .ne. iblock) then
                   print *,'Server ',me,' Error: ',call_marker,' ptr ',
     *                ptr,' memloc ',server_table(c_server_memloc,ptr),
     *                ' should be ',iblock
                   call server_abort_job(server_table, nserver_table)   
               endif
            endif
         else
            if (server_table(c_server_busy_node,ptr) .ne. 0) then
               print *,'Server ',me,' Error: ',call_marker,' ptr ',
     *               ptr,' is not busy, but busy_node ',
     *           ' is set to ',server_table(c_server_busy_node,ptr)
               call server_abort_job(server_table,nserver_table)
            endif
         endif
      enddo

      return
      do i = 1, nserver_table
         memloc = server_table(c_server_memloc,i)
         if (memloc .gt. 0) then
            if (server_table_ptr(memloc) .ne. i) then
               print *,'Server ',me,' table inconsistency at ',
     *           call_marker,' memloc, ptr = ',memloc,i,
     *           ' server_table_ptr = ',server_table_ptr(memloc)
               call server_abort_job(server_table,nserver_table)
            endif

            do j = 1, i-1
               if (server_table(c_server_memloc,j) .eq. memloc) then
                 print *,'Server ',me,' table inconsistency at ',
     *           call_marker,' memloc, ptr = ',memloc,i,
     *           ' server_table_ptr = ',server_table_ptr(memloc),
     *           ' j, server_table(c_server_memloc,j) ',j,
     *            server_table(c_server_memloc,j)

                 print *,'Server ',me,' Dump of ptr ',i,' data ',
     *               (server_table(k,i),k=1,lserver_table_entry)
                 print *,'Server ',me,' Dump of ptr ',j,' data ',
     *               (server_table(k,j),k=1,lserver_table_entry)
                 call server_abort_job(server_table,nserver_table)
               endif 
            enddo
         endif
      enddo

      return
      end
