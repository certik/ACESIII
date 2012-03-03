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
      subroutine process_server_blocks_to_list_msg(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_blocks_to_list message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'block_index.h'
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
      integer i, j, k, istart, nsearch, iblock
      integer*8 indblk, get_index_from_base
      integer*8 ind, imsg
      integer msgbuffer, state, ierr
      integer memloc, diskloc, size, ifile, nxt
      integer request, status(MPI_STATUS_SIZE)
      integer array, nindex, tag, f_form_msg_tag 
      integer msglen
      integer msg(100)
      integer*8 offset, index_file_loc
      integer fh, fhndx, mode
      integer block_index_entry(lblk_index_entry)
      logical flag
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
      save fh, fhndx

#ifdef ALTIX
      dptr = dshptr
#endif

      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then
         array = server_msg(c_msg_array,node)
         if (array .eq. 0) then

c--------------------------------------------------------------------------
c   We are being signalled to open the BLOCKDATA files.  More messages 
c   with array info are coming later.
c--------------------------------------------------------------------------

            mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
            call mpi_file_open(mpi_comm_world, 'BLOCKDATA2', mode,
     *                   MPI_INFO_NULL, fh, ierr)
            offset = 0
            call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION,
     *                       MPI_DOUBLE_PRECISION, 'native',
     *                       MPI_INFO_NULL, ierr)

            call mpi_file_open(mpi_comm_world, 'BLOCK_INDEX2', mode,

     *                   MPI_INFO_NULL, fhndx, ierr)
            offset = 0
            call mpi_file_set_view(fhndx, offset, MPI_INTEGER,
     *                       MPI_INTEGER, 'native',
     *                       MPI_INFO_NULL, ierr)

            server_msg(c_msg_state,node) = null_state
            return 
         else if (array .eq. blocks_list_done) then
        
c----------------------------------------------------------------------------
c   Signal to close the files.
c----------------------------------------------------------------------------

            call mpi_file_sync(fh, ierr)
            call mpi_file_close(fh, ierr)
            call mpi_file_sync(fhndx, ierr)
            call mpi_file_close(fhndx, ierr)
 
            server_msg(c_msg_state,node) = null_state
            return
         endif

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
         index_file_loc = blk_to_list_offset(1,node) 
         offset         = blk_to_list_offset(2,node)

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

               size = server_table(c_server_size,iblock)
               if (server_table(c_server_memloc,iblock) .ne. 0) then
                  
c---------------------------------------------------------------------------
c   Data is found in memory.
c---------------------------------------------------------------------------
   
                  j = server_table(c_server_memloc,iblock)
                  ind = indblk + (j-1) * server_mem_blocksize
               else if (server_table(c_server_diskloc,iblock) .gt. 
     *                      0) then 

c---------------------------------------------------------------------------
c   Read the data from disk into the node's message buffer.
c---------------------------------------------------------------------------

                  diskloc = server_table(c_server_diskloc,iblock)
                  ifile   = server_table(c_server_file, iblock)
                  call f_restoreram(server_unit(ifile), diskloc,
     *                           server_blocksizes(ifile),
     *                           x(imsg), size)
                  ind = imsg
               else
                  print *,'Server ',me,': Error in write_blocks_to_list'
                  print *,'   Array ',array,
     *                    ' erroneous server_table entry'
                  print *,(server_table(i,iblock),i=1,
     *                      lserver_table_entry)
                  call server_abort_job(server_table, nserver_table)
               endif

c-------------------------------------------------------------------------
c   Write block to disk.
c-------------------------------------------------------------------------

               call mpi_file_write_at(fh, offset, x(ind),
     *                     size, MPI_DOUBLE_PRECISION, status,
     *                     ierr)
c               print *,'B2L: offset, size ',offset,size,' xxx ',
c     *                block_energy(x(ind), size)

c--------------------------------------------------------------------------
c   Write the index file data to disk.
c--------------------------------------------------------------------------

               block_index_entry(c_blk_index_array) = array
               block_index_entry(c_blk_index_blkno) = 
     *                      server_table(c_server_iblk,iblock)
               block_index_entry(c_blk_index_size)  = size
               nindex = server_table(c_server_nind,iblock)
               block_index_entry(c_blk_index_nind)  = nindex
               do k = 1, nindex
                  block_index_entry(c_blk_index_bsegs+k-1) = 
     *              server_table(c_server_bsegs+k-1,iblock)
                  block_index_entry(c_blk_index_esegs+k-1) = 
     *              server_table(c_server_esegs+k-1,iblock)
               enddo
 
               call mpi_file_write_at(fhndx, index_file_loc,
     *                     block_index_entry, lblk_index_entry,
     *                     MPI_INTEGER, status,
     *                     ierr)
               index_file_loc = index_file_loc + lblk_index_entry
               offset = offset + size
            endif
         enddo

c--------------------------------------------------------------------------
c   Complete the processing of this message by entering the null state.
c--------------------------------------------------------------------------

1000     continue   
         server_msg(c_msg_state,node) = null_state   ! done
      endif   ! state
      return
      end
