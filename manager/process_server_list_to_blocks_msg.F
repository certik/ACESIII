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
      subroutine process_server_list_to_blocks_msg(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_list_to_blocks message.
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
      integer i, j, istart, nsearch, iblock
      integer*8 indblk, get_index_from_base
      integer*8 ind, imsg
      integer msgbuffer, state, ierr
      integer memloc, diskloc, size, ifile, nxt
      integer status(MPI_STATUS_SIZE)
      integer array, nindex, tag, f_form_msg_tag 
      integer msglen
      integer msg(100)
      integer*8 offset, index_file_loc
      integer fh, fhndx, mode
      integer block_index_entry(lblk_index_entry)
      integer request
      integer ptr, find_server_table_ptr
      integer nind, iblk, nb_total
      integer istat, find_free_diskloc, f_backupram
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

            mode = MPI_MODE_RDONLY
            call mpi_file_open(mpi_comm_world, 'BLOCKDATA', mode,
     *                   MPI_INFO_NULL, fh, ierr)
            offset = 0
            call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION,
     *                       MPI_DOUBLE_PRECISION, 'native',
     *                       MPI_INFO_NULL, ierr)

            call mpi_file_open(mpi_comm_world, 'BLOCK_INDEX', mode,
     *                   MPI_INFO_NULL, fhndx, ierr)
            offset = 0
            call mpi_file_set_view(fhndx, offset, MPI_INTEGER,
     *                       MPI_INTEGER, 'native',
     *                       MPI_INFO_NULL, ierr)

            server_msg(c_msg_state,node) = null_state
            return 
         else if (array .eq. list_to_blocks_done) then
        
c----------------------------------------------------------------------------
c   Signal to close the files.
c----------------------------------------------------------------------------

            call mpi_file_close(fh, ierr)
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
         nb_total       = blk_to_list_offset(3,node)

c---------------------------------------------------------------------------
c   Calculate the index of the message buffer.
c--------------------------------------------------------------------------

         msgbuffer = server_msg(c_msg_msgbuffer,node)
         imsg = indblk + (msgbuffer-1)* server_mem_blocksize

         do iblock = 1, nb_total

c--------------------------------------------------------------------------
c   Read data for this block from the index file.
c--------------------------------------------------------------------------

            call mpi_file_read_at(fhndx, index_file_loc,
     *                     block_index_entry, lblk_index_entry,
     *                     MPI_INTEGER, status,
     *                     ierr)

            index_file_loc = index_file_loc + lblk_index_entry

c----------------------------------------------------------------------------
c   Create dummy info in the message node for use in the table search routine.
c----------------------------------------------------------------------------

            server_msg(c_msg_type,node) =  
     *                                  server_list_to_blocks_msg
            server_msg(c_msg_array,node) =  array   ! use current array
            server_msg(c_msg_source,node) = 0
            server_msg(c_msg_nind,node)   = 
     *                           block_index_entry(c_blk_index_nind)
            nind = server_msg(c_msg_nind,node)
            server_msg(c_msg_iblk,node)   = 
     *                        block_index_entry(c_blk_index_blkno)
            server_msg(c_msg_size,node)   = 
     *                           block_index_entry(c_blk_index_size)
            server_msg(c_msg_state,node)  = begin_state
            server_msg(c_msg_cause,node)  = null_cause
            server_msg(c_msg_memptr,node) = 0
            server_msg(c_msg_flag,node)   = 0
            do i = 1, mx_array_index
               if (i .gt. nind) then
                  server_msg(c_msg_bsegs+i-1,node) = 0
                  server_msg(c_msg_esegs+i-1,node) = 0
               else
                  server_msg(c_msg_bsegs+i-1,node) = 
     *               block_index_entry(c_blk_index_bsegs+i-1)
                  server_msg(c_msg_esegs+i-1,node) =
     *               block_index_entry(c_blk_index_esegs+i-1)
               endif
            enddo

c----------------------------------------------------------------------------
c   Find a match for this data in the server_table.
c----------------------------------------------------------------------------

            ptr = find_server_table_ptr(node,server_table,
     *                                  nserver_table, .false.)
            size = server_msg(c_msg_size,node)

            if (ptr .gt. 0) then
               server_msg(c_msg_stptr,node) = ptr
               server_table(c_server_iblk,ptr) = 
     *                      server_msg(c_msg_iblk,node)

c-------------------------------------------------------------------------
c   Read block from disk.
c-------------------------------------------------------------------------

               call mpi_file_read_at(fh, offset, x(imsg),
     *                  size, MPI_DOUBLE_PRECISION, status,
     *                  ierr)

c---------------------------------------------------------------------------
c   Write the data to the server's SCR file.
c---------------------------------------------------------------------------

               diskloc = server_table(c_server_diskloc,ptr)
               ifile   = server_table(c_server_file,ptr)
               if (diskloc .eq. 0) then
                  diskloc = find_free_diskloc(ifile, server_table,
     *                                   nserver_table)
                  if (diskloc .lt. 0) then
                     next_server_diskloc(ifile) =
     *                 next_server_diskloc(ifile) + 1
                     diskloc = next_server_diskloc(ifile)
                  endif

                  server_table(c_server_diskloc,ptr) = diskloc
               endif

               istat =  f_backupram(server_unit(ifile), diskloc,
     *                           server_blocksizes(ifile),
     *                           x(imsg), size)
               if (istat .ne. 0)
     *            call server_abort_job(server_table,
     *                                  nserver_table)
            endif   ! ptr > 0

            offset = offset + size
         enddo   ! iblock
 
c--------------------------------------------------------------------------
c   Complete the processing of this message by entering the null state.
c--------------------------------------------------------------------------

1000     continue   
         server_msg(c_msg_state,node) = null_state   ! done
      endif   ! state
      return
      end
