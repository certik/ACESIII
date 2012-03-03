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
      subroutine check_for_server_array(node, server_table, 
     *                                    nserver_table)
c---------------------------------------------------------------------------
c   This subroutine performs error checking on the array used in the 
c   server message in "node".
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_stat.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      integer node
      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)

      integer msg_type, ierr
      integer i, array, ix, ptr, ptr2

      if (node .lt. 1 .or. node .gt. nmessage_buffers) then
         print *,'Task ',me,' NODE RANGE ERR: node ',node
         call server_abort_job(server_table, nserver_table)
      endif

      msg_type = server_msg(c_msg_type,node)

      if (msg_type .eq. server_request_msgtype) then

c-------------------------------------------------------------------------
c   REQUEST MESSAGE
c-------------------------------------------------------------------------

         array = server_msg(c_msg_array,node)
         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               ix = i
               go to 100 
            endif
         enddo
         print *,'Server ',me,' REQUEST Cannot find array ',array
         call server_abort_job(server_table, nserver_table)

  100 continue
         if (served_array_status(ix) .eq. 0) then
            served_array_status(ix) = readonly_flag
         else if (served_array_status(ix) .ne. readonly_flag) then
            print *,'Server ',me,' Error: Array ',array,' is ',
     *          ' writeonly, but   REQUEST has been received.'
            print *,'ix, served_array_status(ix) = ',
     *           ix, served_array_status(ix)
            print *,'Server ',me,' Processing node ',node
            call server_abort_job(server_table, nserver_table)
         endif

      else  if (msg_type .eq. server_prepare_msgtype) then

c-------------------------------------------------------------------------
c   PREPARE MESSAGE
c--------------------------------------------------------------------------

         array = server_msg(c_msg_array,node)
         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               ix = i
               go to 200 
            endif
         enddo
         print *,'Server ',me,' PREPARE Cannot find array ',array
         call server_abort_job(server_table, nserver_table)

  200 continue
         if (served_array_status(ix) .eq. 0) then
            served_array_status(ix) = writeonly_flag
         else if (served_array_status(ix) .ne. writeonly_flag) then
            print *,'Server ',me,' Error: Array ',array,' is ',
     *          ' readonly, but a PREPARE has been received.'
            print *,'ix, served_array_status(ix) = ',
     *           ix, served_array_status(ix)
            print *,'Server ',me,' Processing node ',node
            call server_abort_job(server_table, nserver_table)
         endif

      else if (msg_type .eq. server_prepare_increment) then

c--------------------------------------------------------------------------
c   PREPARESUM MESSAGE
c--------------------------------------------------------------------------

         array = server_msg(c_msg_array,node)
         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               ix = i
               go to 300 
            endif
         enddo
         print *,'Server ',me,' PREPARESUM Cannot find array ',array
         call server_abort_job(server_table, nserver_table)

  300 continue
         if (served_array_status(ix) .eq. 0) then
            served_array_status(ix) = writeonly_flag
         else if (served_array_status(ix) .ne. writeonly_flag) then
            print *,'Server ',me,' Error: Array ',array,' is ',
     *          ' readonly, but a PREPARESUM has been received.'
            print *,'Server ',me,' Processing node ',node
            call server_abort_job(server_table, nserver_table)
         endif

      else if (msg_type .eq. server_prequest_msg) then

c-------------------------------------------------------------------------
c   PREQUEST MESSAGE
c-------------------------------------------------------------------------

         array = server_msg(c_msg_array,node)
         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               ix = i
               go to 400 
            endif
         enddo
         print *,'Server ',me,' PREQUEST Cannot find array ',array
         call server_abort_job(server_table, nserver_table)

  400 continue
         if (served_array_status(ix) .eq. 0) then
            served_array_status(ix) = readonly_flag
         else if (served_array_status(ix) .ne. readonly_flag) then
            print *,'Server ',me,' Error: Array ',array,' is ',
     *          ' writeonly, but   PREQUEST has been received.'
            print *,'ix, served_array_status(ix) = ',
     *           ix, served_array_status(ix)
            print *,'Server ',me,' Processing node ',node
            call server_abort_job(server_table, nserver_table)
         endif
      endif

      return
      end
