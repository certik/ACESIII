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
      subroutine server_abort_job(server_table, nserver_table)
c---------------------------------------------------------------------------
c   Abort job from the server side.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)

      integer i, j, k, n, ierr
      integer ptr 

c---------------------------------------------------------------------------
c   Dump the server_table
c---------------------------------------------------------------------------

      print *,'ABORT JOB on server ',me,' SERVER TABLE:'

      print *,'CLEAN BLOCK LIST: clean_block_ptr ',clean_block_ptr
      do i = 1, clean_block_ptr
         j = clean_blocks(i)
         ptr    = server_table_ptr(j) 
         print *,'clean_blocks(i) ',j,' ptr ',ptr,
     *        ' flags ',server_table(c_server_flags,ptr)
      enddo

      do i = 1, nserver_table
         print *,'Server ',me,' entry ',i,': ',(server_table(j,i),
     *        j = 1, lserver_table_entry)
      enddo

      print *,'Server ',me,' served_array_table: '
      print 101,(i, served_array_table(i),
     *           served_array_status(i),i=1,nserved_arrays)
  101 format('   Entry ',i4,' array ',i4,' status ',i8)

      print *,'Server ',me,' Dump of server_table_ptr: '
      do i = 1, (nserver_memblocks + 9)/10
         j = (i-1)*10 + 1
         n = min(nserver_memblocks, j+9)
         print *,'Server ',me,' pointers ',j,' - ',n,': ',
     *           (server_table_ptr(j+k-1),k=1,n-j+1)
      enddo 

c---------------------------------------------------------------------------
c   Dump the message nodes.
c---------------------------------------------------------------------------

      print *,'Server ',me,' Message nodes: head, tail = ',
     *   server_work_list_head, server_work_list_tail
      do i = 1, nmessage_buffers
         print *,'Server ',me,' message node ',i,':',
     *     (server_msg(j,i),j=1,lserver_msg_entry)
      enddo

      call c_flush_stdout()

      call mpi_abort(mpi_comm_world, 10101,ierr)
      return
      end
