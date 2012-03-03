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
      subroutine process_server_blocks_to_list_no_mpi_io(node, 
     *                        server_table, nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_blocks_to_list message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer node
      integer i, j, istart, nsearch, iblock
      integer*8 indblk, get_index_from_base
      integer*8 ind, imsg
      integer msgbuffer, state, ierr
      integer memloc, diskloc, size, ifile, nxt
      integer request, status(MPI_STATUS_SIZE)
      integer array, nindex, tag, f_form_msg_tag 
      integer msglen
      integer msg(100)
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif
      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then
         array = server_msg(c_msg_array,node)

c---------------------------------------------------------------------------
c   Locate the server_table entries for this array.
c---------------------------------------------------------------------------

         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               istart = served_array_entry(i)
               nsearch = served_numblocks(i)
               go to 50
            endif
         enddo

c---------------------------------------------------------------------------
c   Array not found on this server.   We have nothing to do.
c---------------------------------------------------------------------------

         go to 1000
 
   50    continue
         indblk = get_index_from_base(base_mem_addr, x, 2)

c---------------------------------------------------------------------------
c   Calculate the index of the message buffer.
c--------------------------------------------------------------------------

         msgbuffer = server_msg(c_msg_msgbuffer,node)
         imsg = indblk + (msgbuffer-1)* server_mem_blocksize

c---------------------------------------------------------------------------
c   Process each block of the array found in the server_table.
c---------------------------------------------------------------------------

         do iblock = istart, istart+nsearch-1
            if (server_table(c_server_array,iblock) .eq. array) then
               if (server_table(c_server_iblk,iblock) .eq. 
     *                               99999999) then
                  print *,'No message for iblock ',iblock
                  print *,(server_table(i,iblock),
     *                     i=1,lserver_table_entry)
                  call server_abort_job(server_table, nserver_table)
               endif

c---------------------------------------------------------------------------
c   Form a descriptor message and send it to the master.
c---------------------------------------------------------------------------

               msg(1) = server_blocks_to_list_msg 
               msg(2) = array
               nindex = server_table(c_server_nind,iblock)
               msg(3) = nindex
               tag    = f_form_msg_tag()
               msg(4) = tag
               msg(5) = 0
               msg(7) = server_table(c_server_iblk,iblock)
               nxt = 7 + mx_array_index
               do i = 1, mx_array_index
                  if (i .le. nindex) then
                     msg(nxt) = server_table(c_server_bsegs+i-1,iblock)
                     msg(nxt+1)=server_table(c_server_esegs+i-1,iblock)
                  else
                     msg(nxt) = 0
                     msg(nxt+1) = 0
                  endif
                  nxt = nxt + 2
               enddo
               msglen = nxt - 1


               if (server_table(c_server_memloc,iblock) .ne. 0) then
                  
c---------------------------------------------------------------------------
c   Send the data from the memory block to the master.
c---------------------------------------------------------------------------
   
                  j = server_table(c_server_memloc,iblock)
                  ind = indblk + (j-1) * server_mem_blocksize
                  call mpi_isend(x(ind), 
     *                          server_table(c_server_size,iblock),
     *                          mpi_double_precision, 0,
     *                          tag, mpi_comm_world,
     *                          request, ierr)

               else if (server_table(c_server_diskloc,iblock) .gt. 
     *                      0) then 

c---------------------------------------------------------------------------
c   Read the data from disk into the node's message buffer.
c---------------------------------------------------------------------------

                  diskloc = server_table(c_server_diskloc,iblock)
                  size    = server_table(c_server_size,iblock)
                  ifile   = server_table(c_server_file, iblock)
                  call f_restoreram(server_unit(ifile), diskloc,
     *                           server_blocksizes(ifile),
     *                           x(imsg), size)

                  call mpi_isend(x(imsg), size,
     *                          mpi_double_precision, 0,
     *                          tag, mpi_comm_world,
     *                          request, ierr)
               else
                  print *,'Server ',me,': Error in write_blocks_to_list'
                  print *,'   Array ',array,
     *                    ' erroneous server_table entry'
                  print *,(server_table(i,iblock),i=1,
     *                      lserver_table_entry)
                  call server_abort_job(server_table, nserver_table)
               endif
            endif

c-------------------------------------------------------------------------
c   Send the descriptor message.
c-------------------------------------------------------------------------

            call mpi_send(msg, msglen, mpi_integer, 0,
     *                        server_blocks_to_list_msg, 
     *                        mpi_comm_world,
     *                        status, ierr)

c--------------------------------------------------------------------------
c   Wait for the data send to complete.
c--------------------------------------------------------------------------

            call mpi_wait(request, status, ierr)
         enddo
 
c--------------------------------------------------------------------------
c   Complete the processing of this message by entering the null state.
c--------------------------------------------------------------------------

1000     continue   
         server_msg(c_msg_state,node) = null_state   ! done
      endif
      return
      end
