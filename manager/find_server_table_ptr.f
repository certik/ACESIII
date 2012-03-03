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
      integer function find_server_table_ptr(node, server_table,
     *                                       nserver_table, abort_flag)
c----------------------------------------------------------------------------
c   Finds the match for "node" in the server table and returns its
c   pointer.  If no match is found, and the abort_flag = .true.,
c   it is a fatal error and the job is aborted.
c----------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'
      integer node, nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      logical abort_flag

      integer i, j, nind, array, istart, nsearch
      integer msgtype

      array = server_msg(c_msg_array,node)
      nind  = server_msg(c_msg_nind,node) 
      find_server_table_ptr = 0
      if (nind .lt. 1 .or. nind .gt. mx_array_index) then
         print *,'Task ',me,' Error in find_server_table_ptr: nind = ',
     *       nind
         call server_abort_job(server_table, nserver_table)
      endif

      do i = 1, nserved_arrays
         if (served_array_table(i) .eq. array) then
            istart = served_array_entry(i)
            nsearch = served_numblocks(i)
            go to 50
         endif
      enddo

      print *,'Error: Cannot find array ',array,
     *       ' in served array table'
      call server_abort_job(server_table, nserver_table)
   50 continue

      if (istart .lt. 1 .or. istart .gt. nserver_table_entries) then
         print *,'Task ',me,' Error in find_server_table_ptr'
         print *,'istart = ',istart,' nserver_table_entries ',
     *           nserver_table_entries
         call server_abort_job(server_table, nserver_table)
      endif

      if (istart+nsearch-1 .lt. 1 .or. istart+nsearch-1 .gt.
     *          nserver_table_entries) then
         print *,'Task ',me,' Error in find_server_table_ptr'
         print *,'istart = ',istart,' nserver_table_entries ',
     *           nserver_table_entries
         print *,'nsearch = ',nsearch,' istart+nsearch-1 ',
     *           istart+nsearch-1
         call server_abort_job(server_table, nserver_table)
      endif

      msgtype = server_msg(c_msg_type, node)

      do 100 i = istart, istart+nsearch-1
         if (msgtype .eq. server_prequest_msg) then

c--------------------------------------------------------------------------
c   For a partial request, the secondary set of segments is searched.
c--------------------------------------------------------------------------

            do j = 1, nind
               if (server_msg(c_msg_bsegs2+j-1,node) .ne.
     *             server_table(c_server_bsegs+j-1,i)) go to 100
               if (server_msg(c_msg_esegs2+j-1,node) .ne.
     *             server_table(c_server_esegs+j-1,i)) go to 100
            enddo
         else

c---------------------------------------------------------------------------
c   Normal search.
c---------------------------------------------------------------------------

            do j = 1, nind
               if (server_msg(c_msg_bsegs+j-1,node) .ne. 
     *             server_table(c_server_bsegs+j-1,i)) go to 100
               if (server_msg(c_msg_esegs+j-1,node) .ne.
     *             server_table(c_server_esegs+j-1,i)) go to 100
            enddo 
         endif

         find_server_table_ptr = i

         if (array .ne. server_table(c_server_array,i)) then
            print *,'Error: server_table array value doesnt match',
     *             ' expected value.'
            print *,'Actual value = ',server_table(c_server_array,i)
            print *,'Expected value = ',array
            call server_abort_job(server_table, nserver_table)
         endif

         if (msgtype .ne. server_prequest_msg .and.
     *       server_table(c_server_size,i) .ne. 
     *                  server_msg(c_msg_size,node)) then
            print *,'Server ',me,' Size mismatch: ptr ',i,
     *       ' size ',server_table(c_server_size,i),' node ',
     *       node, ' msg_size ',server_msg(c_msg_size,node)
            call server_abort_job(server_table, nserver_table)
         endif
         return
  100 continue

      if (abort_flag) then  
         print *,'Error: Entry for node not found in server_table'
         if (msgtype .eq. server_prequest_msg) then
            print *,'node ',node,' array ',array,' segs ',
     *        (server_msg(c_msg_bsegs2+j-1,node),
     *        server_msg(c_msg_esegs2+j-1,node),j=1,nind)
         else
            print *,'node ',node,' array ',array,' segs ',
     *        (server_msg(c_msg_bsegs+j-1,node),
     *        server_msg(c_msg_esegs+j-1,node),j=1,nind)
         endif
         call server_abort_job(server_table, nserver_table)
      endif   ! abort_flag
      return
      end
